#!/usr/bin/env python3

import numpy as np
import os
import ctypes
from pathlib import Path

class PauliSolver:
    """Python wrapper for C++ PauliSolver class"""
    
    def __init__(self, verbosity=0):
        self.verbosity = verbosity
        self.solver = None
        self.lib = self._compile_and_load()
        self._setup_function_signatures()
    
    def _compile_and_load(self):
        """Compile and load the C++ library"""
        cpp_dir = Path(os.path.dirname(__file__)) / "cpp"
        
        # Compile the wrapper
        os.system(f"g++ -O3 -fPIC -shared -o {cpp_dir}/libpauli_solver.so "
                 f"{cpp_dir}/pauli_solver_wrapper.cpp")
        
        return ctypes.CDLL(f"{cpp_dir}/libpauli_solver.so")
    
    def _setup_function_signatures(self):
        """Set up the C++ function signatures"""
        # Create solver
        self.lib.create_pauli_solver.argtypes = [
            ctypes.c_int,      # nstates
            ctypes.c_int,      # nsingle
            ctypes.POINTER(ctypes.c_double),  # energies
            ctypes.c_double,   # eps1
            ctypes.c_double,   # eps2
            ctypes.c_double,   # eps3
            ctypes.c_double,   # t
            ctypes.c_double,   # W
            ctypes.c_double,   # VBias
            ctypes.c_double,   # GammaS
            ctypes.c_double,   # GammaT
            ctypes.c_double,   # coeffE
            ctypes.c_double,   # coeffT
            ctypes.POINTER(ctypes.c_double),  # lead_mu
            ctypes.POINTER(ctypes.c_double),  # lead_temp
            ctypes.c_int,      # verbosity
        ]
        self.lib.create_pauli_solver.restype = ctypes.c_void_p
        
        # Delete solver
        self.lib.delete_pauli_solver.argtypes = [ctypes.c_void_p]
        
        # Solve master equation
        self.lib.solve_master_equation.argtypes = [ctypes.c_void_p]
        
        # Get results
        self.lib.get_probabilities.argtypes = [ctypes.c_void_p]
        self.lib.get_probabilities.restype = ctypes.POINTER(ctypes.c_double)
        
        self.lib.get_kernel.argtypes = [ctypes.c_void_p]
        self.lib.get_kernel.restype = ctypes.POINTER(ctypes.c_double)
        
        self.lib.get_rhs.argtypes = [ctypes.c_void_p]
        self.lib.get_rhs.restype = ctypes.POINTER(ctypes.c_double)
        
        self.lib.get_pauli_factors.argtypes = [ctypes.c_void_p]
        self.lib.get_pauli_factors.restype = ctypes.POINTER(ctypes.c_double)
    
    def solve(self, energies, nsingle, eps1, eps2, eps3, t, W, VBias,
             GammaS, GammaT, coeffE, coeffT, lead_mu, lead_temp):
        """Solve the master equation using the C++ solver"""
        # Convert numpy arrays to C arrays
        nstates = len(energies)
        
        energies_arr = energies.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lead_mu_arr = lead_mu.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lead_temp_arr = lead_temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        # Create solver instance
        self.solver = self.lib.create_pauli_solver(
            nstates, nsingle,
            energies_arr,
            eps1, eps2, eps3,
            t, W, VBias,
            GammaS, GammaT,
            coeffE, coeffT,
            lead_mu_arr, lead_temp_arr,
            self.verbosity
        )
        
        # Solve master equation
        self.lib.solve_master_equation(self.solver)
        
        # Get results
        probabilities = np.ctypeslib.as_array(
            self.lib.get_probabilities(self.solver),
            shape=(nstates,)
        )
        
        kernel = np.ctypeslib.as_array(
            self.lib.get_kernel(self.solver),
            shape=(nstates, nstates)
        )
        
        rhs = np.ctypeslib.as_array(
            self.lib.get_rhs(self.solver),
            shape=(nstates,)
        )
        
        pauli_factors = np.ctypeslib.as_array(
            self.lib.get_pauli_factors(self.solver),
            shape=(2, nstates, nstates, 2)  # [nleads, nstates, nstates, 2]
        )
        
        return probabilities, kernel, rhs, pauli_factors
    
    def __del__(self):
        """Clean up C++ solver"""
        if self.solver is not None:
            self.lib.delete_pauli_solver(self.solver)

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
