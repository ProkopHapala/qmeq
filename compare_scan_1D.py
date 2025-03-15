#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from sys import path
path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')

import qmeq
from pauli_solver_lib import PauliSolver, calculate_state_energy, calculate_tunneling_amplitudes

# Constants
NSingle = 3  # number of impurity states
NLeads = 2   # number of leads

# Parameters (in meV)
t = 0.0      # direct hopping
W = 20.0     # inter-site coupling
VBias = 0.1  # bias voltage

# Lead parameters
muS = 0.0    # substrate chemical potential
muT = 0.0    # tip chemical potential
Temp = 0.224 # temperature in meV
DBand = 1000.0 # lead bandwidth
GammaS = 0.20 # coupling to substrate
GammaT = 0.05 # coupling to tip

# Tunneling amplitudes
VS = np.sqrt(GammaS/np.pi)  # substrate
VT = np.sqrt(GammaT/np.pi)  # tip

# Position-dependent coefficients
coeffE = 0.4
coeffT = 0.3

# Lead couplings
TLeads = {(0,0): VS, (0,1): VS, (0,2): VS,
          (1,0): VT, (1,1): coeffT*VT, (1,2): coeffT*VT}

def initialize_qmeq_solver():
    """Initialize QmeQ solver once"""
    print('\n### Initializing QmeQ Pauli solver')
    
    hsingle = {(0,0): 0.0, (0,1): t, (0,2): t,
               (1,1): 0.0, (1,2): t,
               (2,2): 0.0}
    coulomb = {(0,1,1,0): W, (1,2,2,1): W, (0,2,2,0): W}
    mu_L = {0: muS, 1: muT + VBias}
    Temp_L = {0: Temp, 1: Temp}

    qmeq_system = qmeq.Builder(NSingle, hsingle, coulomb, NLeads, TLeads, mu_L, Temp_L, DBand,
                             kerntype='Pauli', indexing='Lin', itype=0, symq=True, solmethod='solve', mfreeq=0)
    
    # Set verbosity after creation
    qmeq_system.appr.verbosity = verbosity
    qmeq_system.verbosity = verbosity
    
    print(f'QmeQ solver initialized with verbosity: {verbosity}')
    print(f'QmeQ system params: NSingle={NSingle}, NLeads={NLeads}, DBand={DBand}')
    return qmeq_system

def initialize_cpp_solver():
    """Initialize C++ solver once"""
    pauli = PauliSolver(verbosity=verbosity)
    NStates = 2**NSingle
    
    # Create constant parts of the solver
    tunneling_amplitudes = calculate_tunneling_amplitudes(NLeads, NStates, NSingle, TLeads)
    lead_mu = np.array([muS, muT + VBias])
    lead_temp = np.array([Temp, Temp])
    lead_gamma = np.array([GammaS, GammaT])
    
    # Create solver instance with dummy energies
    solver = pauli.create_solver(NStates, NLeads, np.zeros(NStates), tunneling_amplitudes, lead_mu, lead_temp, lead_gamma, verbosity)
    
    return pauli, solver, NStates

def calculate_energies(eps1, eps2, eps3, NSingle, W, VBias, coeffE, t):
    """Calculate many-body state energies for QmeQ solver"""
    return np.array([calculate_state_energy(i, NSingle, eps1, eps2, eps3, W, VBias=VBias, coeffE=coeffE, t=t) for i in range(2**NSingle)])

def run_qmeq_solver(eps1, eps2, eps3, qmeq_system):
    """Run QmeQ solver with updated Hamiltonian"""
    print(f'\nRunning QmeQ solver for eps=({eps1:.2f}, {eps2:.2f}, {eps3:.2f})')
    
    # Calculate new energies
    energies = calculate_energies(eps1, eps2, eps3, NSingle, W, VBias, coeffE, t)
    
    # Update Hamiltonian
    qmeq_system.Ea = energies
    
    if verbosity > 0:
        np.set_printoptions(linewidth=np.inf, suppress=True)
        print('Energies QmeQ:', np.array2string(qmeq_system.Ea, separator=', ', prefix='Energies QmeQ: '))
    
    qmeq_system.solve()
    
    if verbosity > 0:
        print('Probabilities QmeQ:', np.array2string(qmeq_system.phi0, separator=', ', prefix='Probabilities QmeQ: '))
        print('Kernel QmeQ:')
        print(np.array2string(qmeq_system.kern, separator=', ', prefix='Kernel QmeQ: '))
    
    print('QmeQ current:', qmeq_system.current[1])
    return qmeq_system.current[1]

def run_cpp_solver(eps1, eps2, eps3, pauli, solver, NStates):
    """Run C++ solver with updated energies"""
    print(f'\nRunning C++ solver for eps=({eps1:.2f}, {eps2:.2f}, {eps3:.2f})')
    
    # Calculate new energies
    energies = calculate_energies(eps1, eps2, eps3, NSingle, W, VBias, coeffE, t)
    
    if verbosity > 0:
        np.set_printoptions(linewidth=np.inf, suppress=True)
        print('Energies C++:', np.array2string(energies, separator=', ', prefix='Energies C++: '))
    
    tunneling_amplitudes = calculate_tunneling_amplitudes(NLeads, NStates, NSingle, TLeads)
    lead_mu = np.array([muS, muT + VBias])
    lead_temp = np.array([Temp, Temp])
    lead_gamma = np.array([GammaS, GammaT])
    
    solver = pauli.create_solver(NStates, NLeads, energies, tunneling_amplitudes, lead_mu, lead_temp, lead_gamma, verbosity)
    pauli.solve(solver)
    
    if verbosity > 0:
        kernel = pauli.get_kernel(solver, NStates)
        probabilities = pauli.get_probabilities(solver, NStates)
        print('Probabilities C++:', np.array2string(probabilities, separator=', ', prefix='Probabilities C++: '))
        print('Kernel C++:')
        print(np.array2string(kernel, separator=', ', prefix='Kernel C++: '))
    
    current = pauli.calculate_current(solver, 1)
    pauli.cleanup(solver)
    print('C++ current:', current)
    return current

def scan_QmeQ(eps, bPrint=False):
    if bPrint:
        print("\n####################################################################")
        print("scan_QmeQ")
        print("######################################################################")
    qmeq_system = initialize_qmeq_solver()
    res = []
    for eps1, eps2, eps3 in eps:
        qmeq_current = run_qmeq_solver(eps1, eps2, eps3, qmeq_system)
        res.append(qmeq_current)
        if bPrint:
            print(f"{eps1:.2f} {eps2:.2f} {eps3:.2f} -> {qmeq_current:.2f}")
    return res

def scan_cpp(eps, bPrint=False):
    if bPrint:
        print("\n####################################################################")
        print("scan_cpp")
        print("######################################################################")
    pauli, cpp_solver, NStates = initialize_cpp_solver()
    res = []
    for eps1, eps2, eps3 in eps:
        cpp_current = run_cpp_solver(eps1, eps2, eps3, pauli, cpp_solver, NStates)
        res.append(cpp_current)
        if bPrint:
            print(f"{eps1:.2f} {eps2:.2f} {eps3:.2f} -> {cpp_current:.2f}")
    return res
    
if __name__ == "__main__":
    # Define energy range
    bPrint = True
    nstep = 3
    eps = np.zeros( (nstep,3) )
    ts = np.linspace(0, 1, nstep)
    eps[:,0] = 0.1+ts
    eps[:,1] = 0.2+ts*2.
    eps[:,2] = 0.3+ts*.3

    verbosity = 1
    
    # Run scan
    qmeq_results = scan_QmeQ(eps, bPrint)
    cpp_results  = scan_cpp(eps, bPrint)
            
    plt.figure(figsize=(10, 6))
    plt.plot(ts, qmeq_results, 'b-', label='QmeQ Pauli')
    plt.plot(ts, cpp_results, 'r--', label='C++ Pauli')
    plt.xlabel('Onsite Energy (meV)')
    plt.ylabel('Current (nA)')
    plt.title('Solver Comparison for 1D Energy Scan')
    plt.legend()
    plt.grid(True)
    plt.show()