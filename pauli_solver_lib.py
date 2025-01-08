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
        if self.verbosity > 0:
            print("DEBUG: pauli_solver_lib.py tunneling amplitudes:")
            print(tunneling_amplitudes)
            
        return self.lib.create_pauli_solver(
            nstates, nleads,
            _np_as(energies, c_double_p),
            _np_as(tunneling_amplitudes.ravel(), c_double_p),
            _np_as(lead_mu, c_double_p),
            _np_as(lead_temp, c_double_p),
            _np_as(lead_gamma, c_double_p),
            verbosity
        )
    
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
    diff = state1 ^ state2
    return diff == (1 << site)

def calculate_tunneling_amplitudes(nleads, nstates, nsingle, vs, vt, coeff_t):
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
    
    for lead in range(nleads):
        v = vs if lead == 0 else vt
        
        for i in range(nstates):
            for j in range(nstates):
                # States must differ by exactly one electron
                diff = count_electrons(i) - count_electrons(j)
                if abs(diff) != 1:
                    continue
                
                # Check if transition is valid (only one site changes)
                for site in range(nsingle):
                    if is_valid_transition(i, j, site):
                        # Apply position-dependent coupling for tip
                        if lead == 1:  # Tip
                            coeff = 1.0 if site == 0 else coeff_t
                            tunneling_amplitudes[lead, j, i] = v * coeff
                            print(f"DEBUG: Python l:{lead} i:{i} j:{j} site:{site} v:{v:.6f} coeff:{coeff:.6f} amplitude:{tunneling_amplitudes[lead,j,i]:.6f}")
                        else:
                            tunneling_amplitudes[lead, j, i] = v
                            print(f"DEBUG: Python l:{lead} i:{i} j:{j} site:{site} v:{v:.6f} amplitude:{tunneling_amplitudes[lead,j,i]:.6f}")
                        break
    
    return tunneling_amplitudes
