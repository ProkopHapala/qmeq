#!/usr/bin/env python3

import numpy as np
import os
import ctypes
from cpp_utils_ import compile_lib, work_dir, _np_as, c_double_p, c_int_p

class PauliSolver:
    """Python wrapper for C++ PauliSolver class"""
    
    def __init__(self):
        """Initialize the solver by loading the C++ library"""
        self.lib = self._compile_and_load()
        self._setup_function_signatures()
        
    def _compile_and_load(self):
        """Compile and load the C++ library"""
        cpp_dir = os.path.join(work_dir(__file__), 'cpp')
        compile_lib('pauli_solver_wrapper', FFLAGS="-std=c++17 -O3 -fPIC", LFLAGS="", path=cpp_dir, clean=True)
        return ctypes.CDLL(os.path.join(cpp_dir, 'pauli_solver_wrapper.so'))
        
    def _setup_function_signatures(self):
        """Set up C++ function signatures"""
        self.lib.create_pauli_solver.argtypes = [
            ctypes.c_int, ctypes.c_int,
            c_double_p, c_double_p, c_double_p, c_double_p, c_double_p
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
    
    def create_solver(self, nstates, nleads, energies, tunneling_amplitudes, lead_mu, lead_temp, lead_gamma):
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
        return self.lib.create_pauli_solver(
            nstates, nleads,
            _np_as(energies, c_double_p),
            _np_as(tunneling_amplitudes.ravel(), c_double_p),
            _np_as(lead_mu, c_double_p),
            _np_as(lead_temp, c_double_p),
            _np_as(lead_gamma, c_double_p)
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

def calculate_state_energy(state, nsingle, eps1, eps2, eps3, W):
    """Calculate energy of a given state
    
    Args:
        state (int): State index
        nsingle (int): Number of single-particle states
        eps1, eps2, eps3 (float): Site energies
        W (float): Inter-site coupling
        
    Returns:
        float: State energy
    """
    energy = 0.0
    # Single-particle energies
    for i in range(nsingle):
        if state & (1 << i):
            energy += eps1 if i == 0 else (eps2 if i == 1 else eps3)
    # Inter-site coupling
    for i in range(nsingle-1):
        if (state & (1 << i)) and (state & (1 << (i+1))):
            energy += W
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
                        else:
                            tunneling_amplitudes[lead, j, i] = v
                        break
    
    return tunneling_amplitudes
