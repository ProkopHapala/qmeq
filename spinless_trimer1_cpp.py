#!/usr/bin/env python3

import numpy as np
import os
import ctypes
from cpp_utils_ import compile_lib, work_dir

# Compile the C++ code
def compile_cpp_solver():
    cpp_dir = os.path.join(work_dir(__file__), 'cpp')
    compile_lib('pauli_solver_wrapper',
               FFLAGS="-std=c++17 -O3 -fPIC",
               LFLAGS="",
               path=cpp_dir,
               clean=True)
    return ctypes.CDLL(os.path.join(cpp_dir, 'pauli_solver_wrapper.so'))

# Load the C++ library
lib = compile_cpp_solver()

# Set up function signatures
lib.create_pauli_solver.argtypes = [
    ctypes.c_int, ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64),
    np.ctypeslib.ndpointer(dtype=np.float64)
]
lib.create_pauli_solver.restype = ctypes.c_void_p

lib.solve_pauli.argtypes = [ctypes.c_void_p]
lib.solve_pauli.restype = None

lib.get_kernel.argtypes = [
    ctypes.c_void_p,
    np.ctypeslib.ndpointer(dtype=np.float64)
]
lib.get_kernel.restype = None

lib.get_probabilities.argtypes = [
    ctypes.c_void_p,
    np.ctypeslib.ndpointer(dtype=np.float64)
]
lib.get_probabilities.restype = None

lib.delete_pauli_solver.argtypes = [ctypes.c_void_p]
lib.delete_pauli_solver.restype = None

# System parameters (same as in spinless_trimer1.py)
NSingle = 3  # number of impurity states
NStates = 2**NSingle
NLeads = 2  # number of leads

# Parameters (in meV)
eps1 = eps2 = eps3 = -10.0
t = 0.0    # direct hopping
W = 20.0   # inter-site coupling
VBias = 0.0  # bias voltage

# Lead parameters
muS = 0.0    # substrate chemical potential
muT = 0.0    # tip chemical potential
Temp = 0.224  # temperature in meV
GammaS = 0.20  # coupling to substrate
GammaT = 0.05  # coupling to tip
VS = np.sqrt(GammaS/np.pi)  # tunneling amplitude substrate
VT = np.sqrt(GammaT/np.pi)  # tunneling amplitude tip

# Coefficients that simulate indirectly the position of the tip
coeffE = 0.4
coeffT = 0.3

def count_electrons(state):
    return bin(state).count('1')

def calculate_state_energy(state):
    energy = 0.0
    # Single-particle energies
    for i in range(NSingle):
        if state & (1 << i):
            energy += eps1 if i == 0 else (eps2 if i == 1 else eps3)
    # Inter-site coupling
    for i in range(NSingle-1):
        if (state & (1 << i)) and (state & (1 << (i+1))):
            energy += W
    return energy

def is_valid_transition(state1, state2, site):
    diff = state1 ^ state2
    return diff == (1 << site)

# Calculate energies for each state
energies = np.array([calculate_state_energy(i) for i in range(NStates)])

# Calculate tunneling amplitudes
tunneling_amplitudes = np.zeros((NLeads, NStates, NStates))

# For each lead and each pair of states
for lead in range(NLeads):
    v = VS if lead == 0 else VT
    
    for i in range(NStates):
        for j in range(NStates):
            # States must differ by exactly one electron
            diff = count_electrons(i) - count_electrons(j)
            if abs(diff) != 1:
                continue
            
            # Check if transition is valid (only one site changes)
            for site in range(NSingle):
                if is_valid_transition(i, j, site):
                    # Apply position-dependent coupling for tip
                    if lead == 1:  # Tip
                        coeff = 1.0 if site == 0 else coeffT
                        tunneling_amplitudes[lead, j, i] = v * coeff
                    else:
                        tunneling_amplitudes[lead, j, i] = v
                    break

# Lead parameters
lead_mu = np.array([muS, muT + VBias])
lead_temp = np.array([Temp, Temp])
lead_gamma = np.array([GammaS, GammaT])

# Create and run solver
solver = lib.create_pauli_solver(
    NStates, NLeads,
    energies,
    tunneling_amplitudes.ravel(),
    lead_mu,
    lead_temp,
    lead_gamma
)

# Solve the master equation
lib.solve_pauli(solver)

# Get results
kernel = np.zeros((NStates, NStates))
probabilities = np.zeros(NStates)

lib.get_kernel(solver, kernel)
lib.get_probabilities(solver, probabilities)

# Print results
print("State energies (meV):")
print(energies)

print("\nTunneling amplitudes for lead 0 (substrate):")
print(tunneling_amplitudes[0])
print("\nTunneling amplitudes for lead 1 (tip):")
print(tunneling_amplitudes[1])

print("\nKernel matrix:")
print(kernel)
print("\nProbabilities:")
print(probabilities)

# Cleanup
lib.delete_pauli_solver(solver)
