#!/usr/bin/env python3

import numpy as np
import os
import ctypes
from cpp_utils_ import compile_lib, work_dir, _np_as, c_double_p, c_int_p

class PauliSolver:
    """Python wrapper for C++ PauliSolver class"""
    
    def __init__(self, verbosity=0):
        """Initialize the solver by loading the C++ library"""
        self.lib = self._compile_and_load()
        self._setup_function_signatures()
        self.verbosity = verbosity
        
    def _compile_and_load(self):
        """Compile and load the C++ library"""
        cpp_dir = os.path.join(work_dir(__file__), 'cpp')
        compile_lib('pauli_solver_wrapper', FFLAGS="-std=c++17 -O3 -fPIC", LFLAGS="", path=cpp_dir, clean=True)
        return ctypes.CDLL(os.path.join(cpp_dir, 'pauli_solver_wrapper.so'))
        
    def _setup_function_signatures(self):
        """Set up C++ function signatures"""
        self.lib.create_pauli_solver.argtypes = [
            ctypes.c_int, ctypes.c_int,
            c_double_p, c_double_p, c_double_p, c_double_p, c_double_p, ctypes.c_int
        ]
        self.lib.create_pauli_solver.restype  = ctypes.c_void_p
        
        self.lib.solve_pauli.argtypes = [ctypes.c_void_p]
        self.lib.solve_pauli.restype  = None
        
        self.lib.get_kernel.argtypes = [ctypes.c_void_p, c_double_p]
        self.lib.get_kernel.restype  = None
        
        self.lib.get_probabilities.argtypes = [ctypes.c_void_p, c_double_p]
        self.lib.get_probabilities.restype  = None
        
        self.lib.calculate_current.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self.lib.calculate_current.restype  = ctypes.c_double
        
        self.lib.delete_pauli_solver.argtypes = [ctypes.c_void_p]
        self.lib.delete_pauli_solver.restype  = None
        print("PauliSolver::_setup_function_signatures() DONE")
    
    def create_solver(self, nstates, nleads, energies, tunneling_amplitudes, lead_mu, lead_temp, lead_gamma, verbosity=0):
        """Create a new PauliSolver instance
        
        Args:
            nstates (int): Number of states
            nleads (int): Number of leads
            energies (np.ndarray): State energies
            tunneling_amplitudes (np.ndarray): Tunneling amplitudes (nleads, nstates, nstates)
            lead_mu (np.ndarray): Chemical potentials for each lead
            lead_temp (np.ndarray): Temperatures for each lead
            lead_gamma (np.ndarray): Coupling strengths for each lead
            
        Returns:
            solver: Handle to C++ solver instance
        """
        # Ensure arrays are C-contiguous and in the correct format
        energies = np.ascontiguousarray(energies, dtype=np.float64)
        tunneling_amplitudes = np.ascontiguousarray(tunneling_amplitudes.transpose(0, 2, 1), dtype=np.float64)
        lead_mu    = np.ascontiguousarray(lead_mu,    dtype=np.float64)
        lead_temp  = np.ascontiguousarray(lead_temp,  dtype=np.float64)
        lead_gamma = np.ascontiguousarray(lead_gamma, dtype=np.float64)
        
        #if self.verbosity > 0:
        #    print("\nDEBUG: pauli_solver_lib.py tunneling amplitudes before C++ (after transpose):")
        #    print(tunneling_amplitudes)
        
        # Create solver
        solver = self.lib.create_pauli_solver(
            nstates, nleads,
            _np_as(energies, c_double_p), _np_as(tunneling_amplitudes, c_double_p),
            _np_as(lead_mu, c_double_p), _np_as(lead_temp, c_double_p), _np_as(lead_gamma, c_double_p),
            self.verbosity
        )
        return solver
    
    def solve(self, solver):
        """Solve the master equation"""
        self.lib.solve_pauli(solver)
    
    def get_kernel(self, solver, nstates):
        """Get the kernel matrix"""
        kernel = np.zeros((nstates, nstates))
        self.lib.get_kernel(solver, _np_as(kernel, c_double_p))
        return kernel
    
    def get_probabilities(self, solver, nstates):
        """Get state probabilities"""
        probabilities = np.zeros(nstates)
        self.lib.get_probabilities(solver, _np_as(probabilities, c_double_p))
        return probabilities
    
    def calculate_current(self, solver, lead):
        """Calculate current through a specific lead"""
        return self.lib.calculate_current(solver, lead)
    
    def cleanup(self, solver):
        """Clean up C++ solver instance"""
        self.lib.delete_pauli_solver(solver)

def count_electrons(state):
    """Count number of electrons in a state"""
    return bin(state).count('1')

def calculate_state_energy(state, nsingle, eps1, eps2, eps3, W, VBias=0.0, coeffE=0.0, t=0.0):
    """Calculate energy of a given state
    
    Args:
        state (int): State index
        nsingle (int): Number of single-particle states
        eps1, eps2, eps3 (float): Site energies
        W (float): Inter-site coupling
        VBias (float): Bias voltage
        coeffE (float): Position-dependent coefficient for bias
        t (float): Direct hopping
        
    Returns:
        float: State energy
    """
    energy = 0.0
    
    # Convert state to binary list for easier manipulation
    state_list = [(state >> i) & 1 for i in range(nsingle)]
    
    # Single-particle energies
    for i in range(nsingle):
        if state_list[i]:
            if i == 0:
                energy += eps1 - coeffE*VBias  # Include bias voltage for first site
            elif i == 1:
                energy += eps2
            else:
                energy += eps3
    
    # Direct hopping terms (from single-particle Hamiltonian)
    # Format: {(i,j): t} for hopping between states i and j
    hopping = {(0,1): t, (0,2): t, (1,2): t}
    for (i,j), hop in hopping.items():
        # Only consider transitions where one state is occupied and the other is empty
        if state_list[i] == 1 and state_list[j] == 0:
            # Calculate fermion sign
            fsign = (-1)**(sum(state_list[0:i]) + sum(state_list[0:j])) * (+1 if i > j else -1)
            # Create the new state
            new_state_list = state_list.copy()
            new_state_list[i] = 0
            new_state_list[j] = 1
            # Convert back to integer for comparison
            new_state = sum(b << i for i, b in enumerate(new_state_list))
            if new_state < state:  # Only add term once
                energy += hop * fsign
    
    # Coulomb interaction terms
    # Format: {(i,j,j,i): W} for interaction between sites i and j
    coulomb = {(0,1,1,0): W, (1,2,2,1): W, (0,2,2,0): W}
    for (i,j,_,_), U in coulomb.items():
        if state_list[i] and state_list[j]:
            energy += U
    return energy

def is_valid_transition(state1, state2, site):
    """Check if transition between states is valid
    
    Args:
        state1, state2 (int): State indices
        site (int): Site index
        
    Returns:
        bool: True if transition is valid
    """
    # Convert states to binary
    bin1 = [int(x) for x in format(state1, '03b')]
    bin2 = [int(x) for x in format(state2, '03b')]
    
    # States should differ only at the specified site
    for i in range(len(bin1)):
        if i == site:
            if bin1[i] == bin2[i]:  # Should be different at site
                return False
        else:
            if bin1[i] != bin2[i]:  # Should be same elsewhere
                return False
    return True

def calculate_tunneling_amplitudes(NLeads, NStates, NSingle, TLeads):
    """Calculate tunneling amplitudes between states using TLeads dictionary format
    
    Args:
        NLeads (int): Number of leads
        NStates (int): Number of states 
        NSingle (int): Number of single-particle states
        TLeads (dict): Dictionary with format {(lead,state): amplitude} defining tunneling amplitudes
        
    Returns:
        np.ndarray: Tunneling amplitudes (nleads, nstates, nstates)
    """
    tunneling_amplitudes = np.zeros((NLeads, NStates, NStates))
    
    # First, group states by charge number for debugging
    states_by_charge = [[] for _ in range(NSingle + 1)]
    for i in range(NStates):
        binary = format(i, f'0{NSingle}b')
        charge = sum(int(bit) for bit in binary)
        states_by_charge[charge].append(i)
        
    # Iterate over all many-body states
    for j1 in range(NStates):
        state = [int(x) for x in format(j1, f'0{NSingle}b')]
        
        # Iterate over all single-particle states
        for j2 in range(NSingle):
            # Calculate fermion sign for added/removed electron
            fsign = np.power(-1, sum(state[0:j2]))
            
            # For each lead
            for lead in range(NLeads):
                if (lead, j2) in TLeads:
                    tamp = TLeads[(lead, j2)]
                    
                    if state[j2] == 0:  # Can add an electron
                        statep = list(state)
                        statep[j2] = 1
                        ind = int(''.join(map(str, statep)), 2)
                        tunneling_amplitudes[lead, ind, j1] = fsign * tamp
                    else:  # Can remove an electron
                        statep = list(state)
                        statep[j2] = 0
                        ind = int(''.join(map(str, statep)), 2)
                        tunneling_amplitudes[lead, ind, j1] = fsign * np.conj(tamp)

    #print("\nDEBUG: Final tunneling amplitudes matrix in PauliSolver::calculate_tunneling_amplitudes() of pauli_solver_lib.py :")
    #for lead in range(NLeads):
    #    print(f"\nLead {lead}:")
    #    print(tunneling_amplitudes[lead])

    return tunneling_amplitudes

def calculate_tunneling_amplitudes_old(nleads, nstates, nsingle, vs, vt, coeff_t):
    """Calculate tunneling amplitudes between states
    
    Args:
        nleads (int): Number of leads
        nstates (int): Number of states
        nsingle (int): Number of single-particle states
        vs (float): Substrate tunneling strength
        vt (float): Tip tunneling strength
        coeff_t (float): Position-dependent coupling coefficient for tip
        
    Returns:
        np.ndarray: Tunneling amplitudes (nleads, nstates, nstates)
    """
    tunneling_amplitudes = np.zeros((nleads, nstates, nstates))
    
    print("\nDEBUG: Python tunneling amplitudes calculation:")
    print(f"vs={vs:.6f}, vt={vt:.6f}, coeff_t={coeff_t:.6f}")
    
    # First, group states by charge number
    states_by_charge = [[] for _ in range(nsingle + 1)]
    for i in range(nstates):
        # Get binary representation and count electrons
        binary = format(i, f'0{nsingle}b')
        charge = sum(int(bit) for bit in binary)
        states_by_charge[charge].append(i)
    
    # Print state grouping for debugging
    print("\nDEBUG: States grouped by charge (our ordering):")
    for charge, states in enumerate(states_by_charge):
        state_strings = [format(state, f'0{nsingle}b') for state in states]
        print(f"Charge {charge}: states {states} (binary: {state_strings})")
        
    # Iterate over all many-body states
    for j1 in range(nstates):
        # Get binary representation of state
        state = [int(x) for x in format(j1, f'0{nsingle}b')]
        
        # Iterate over all single-particle states
        for j2 in range(nsingle):
            # Calculate fermion sign for added/removed electron
            fsign = np.power(-1, sum(state[0:j2]))
            
            # For each lead
            for lead in range(nleads):
                v = vs if lead == 0 else vt
                # For tip, coupling depends on which site changes
                if lead == 1:  # Tip
                    # Apply coeff_t to all sites except the last one
                    coeff = coeff_t if j2 < nsingle-1 else 1.0
                    tamp = v * coeff
                else:  # Substrate
                    tamp = v
                
                if state[j2] == 0:  # Can add an electron
                    # Create new state with electron added at j2
                    statep = list(state)
                    statep[j2] = 1
                    # Convert binary state back to index
                    ind = int(''.join(map(str, statep)), 2)
                    # Add tunneling amplitude (final_state, initial_state)
                    tunneling_amplitudes[lead, ind, j1] = fsign * tamp
                else:  # Can remove an electron
                    # Create new state with electron removed at j2
                    statep = list(state)
                    statep[j2] = 0
                    # Convert binary state back to index
                    ind = int(''.join(map(str, statep)), 2)
                    # Add tunneling amplitude (final_state, initial_state)
                    tunneling_amplitudes[lead, ind, j1] = fsign * np.conj(tamp)
    
    #print("\nDEBUG: Final tunneling amplitudes matrix in PauliSolver::calculate_tunneling_amplitudes() of pauli_solver_lib.py :")
    #for lead in range(nleads):
    #    print(f"\nLead {lead}:")
    #    print(tunneling_amplitudes[lead])
    return tunneling_amplitudes
