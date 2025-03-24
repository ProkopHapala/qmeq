#!/usr/bin/env python3

# Set up ASan preloading before any imports
import os
bASAN = True
if bASAN:
    # Get ASan library path
    asan_lib = os.popen('gcc -print-file-name=libasan.so').read().strip()
    print("Preloading ASan library: ", asan_lib)
    # Set LD_PRELOAD environment variable
    os.environ['LD_PRELOAD'] = asan_lib
    os.environ['ASAN_OPTIONS'] = 'detect_leaks=0'

import sys
sys.stdout = sys.stderr = open(sys.stdout.fileno(), mode='w', buffering=1)

import numpy as np
import matplotlib.pyplot as plt
from sys import path
#path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')

import qmeq
from qmeq import config
from qmeq import indexing as qmqsi
from qmeq.config import verb_print_

import pauli_solver_lib as psl
from pauli_solver_lib import PauliSolver

# setup numpy print options to infinite line length
np.set_printoptions(linewidth=256, suppress=True)

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

def build_hamiltonian(eps1, eps2, eps3, t, W):
    verb_print_(1,"\n#### Building Hamiltonian: eps: ", [eps1, eps2, eps3], " t: ", t, " W: ", W)
    # One-particle Hamiltonian
    hsingle = {(0,0): eps1, (0,1): t, (0,2): t,
               (1,1): eps2, (1,2): t,
               (2,2): eps3}
    
    # Two-particle Hamiltonian: inter-site coupling
    coulomb = {(0,1,1,0): W,
               (1,2,2,1): W,
               (0,2,2,0): W}
    
    return hsingle, coulomb

def build_leads(muS, muT, Temp, VS, VT, coeffT, VBias):
    # Leads: substrate (S) and scanning tip (T)
    mu_L   = {0: muS, 1: muT + VBias}
    Temp_L = {0: Temp, 1: Temp}
    # Coupling between leads (1st number) and impurities (2nd number)
    TLeads = {(0,0): VS,         # S <-- 1
              (0,1): VS,         # S <-- 2
              (0,2): VS,         # S <-- 3
              (1,0): VT,         # T <-- 1
              (1,1): coeffT*VT,  # T <-- 2
              (1,2): coeffT*VT}  # T <-- 3
    
    return mu_L, Temp_L, TLeads

def run_QmeQ_solver(eps1, eps2, eps3):
    """Run QmeQ solver with the given parameters"""
    verb_print_(1,  '\n### Running QmeQ Pauli solver')
    
    mu_L, Temp_L, TLeads = build_leads(muS, muT, Temp, VS, VT, coeffT, VBias)
    Hsingle, Hcoulomb    = build_hamiltonian(eps1, eps2, eps3, t, W)
    
    try:
        config.verbosity = verbosity
        system = qmeq.Builder(NSingle, Hsingle, Hcoulomb, NLeads, TLeads, mu_L, Temp_L, DBand,  kerntype='Pauli', indexing='Lin', itype=0, symq=True, solmethod='solve', mfreeq=0)
        system.appr.verbosity = verbosity  # Set verbosity after instance creation
        system.verbosity = verbosity
        system.solve()
    except Exception as e:
        print(f"Error running QmeQ solver: {e}")
        return None
    
    if verbosity > 0:
        print("QmeQ energies:", system.Ea)
        print("QmeQ probabilities:", system.phi0)
        print("QmeQ current:", system.current[1])
    
    return system.current[1]

# def initialize_qmeq_solver():
#     """Initialize QmeQ solver once"""
#     verb_print_(1,  '\n### Initializing QmeQ Pauli solver')
    
#     hsingle = {(0,0): 0.0, (0,1): t, (0,2): t,
#                (1,1): 0.0, (1,2): t,
#                (2,2): 0.0}
#     coulomb = {(0,1,1,0): W, (1,2,2,1): W, (0,2,2,0): W}
#     mu_L = {0: muS, 1: muT + VBias}
#     Temp_L = {0: Temp, 1: Temp}

#     qmeq_system = qmeq.Builder(NSingle, hsingle, coulomb, NLeads, TLeads, mu_L, Temp_L, DBand, kerntype='Pauli', indexing='Lin', itype=0, symq=True, solmethod='solve', mfreeq=0)
    
#     # Set verbosity after creation
#     qmeq_system.appr.verbosity = verbosity
#     qmeq_system.verbosity = verbosity
    
#     print(f'QmeQ solver initialized with verbosity: {verbosity}')
#     print(f'QmeQ system params: NSingle={NSingle}, NLeads={NLeads}, DBand={DBand}')
#     return qmeq_system

def run_cpp_solver(pauli,eps1, eps2, eps3):
    """Run C++ solver with the given parameters"""
    verb_print_(1,  '\n### Running C++ Pauli solver')
    
    NStates = 2**NSingle
    
    mu_L, Temp_L, TLeads = build_leads(muS, muT, Temp, VS, VT, coeffT, VBias)
    Hsingle, Hcoulomb    = build_hamiltonian(eps1, eps2, eps3, t, W)
    
    lead_mu = np.array([muS, muT + VBias])
    lead_temp = np.array([Temp, Temp])
    lead_gamma = np.array([GammaS, GammaT])
    
    # Convert TLeads dictionary to matrix form
    TLeads_ = np.zeros((NLeads, NSingle))
    for k, v in TLeads.items():
        TLeads_[k[0], k[1]] = v
    
    # Convert Hsingle dictionary to matrix form
    Hsingle_ = np.zeros((NSingle, NSingle))
    for k, v in Hsingle.items():
        Hsingle_[k[0], k[1]] = v
    
    if verbosity > 0:
        print("\nHsingle:"); print(Hsingle_)
        print("\nTLeads:"); print(TLeads_)
    
    # State ordering
    state_order = [0, 4, 2, 6, 1, 5, 3, 7]
    state_order = np.array(state_order, dtype=np.int32)
    
    # Create and run solver
    solver = pauli.create_pauli_solver_new(NStates, NLeads, Hsingle_, W, TLeads_, lead_mu, lead_temp, lead_gamma, state_order, verbosity)
    
    # Get energies before solving
    energies = pauli.get_energies(solver, NStates)
    if verbosity > 0:
        print("C++ energies:", energies)
    
    pauli.solve(solver)
    
    if verbosity > 0:
        kernel = pauli.get_kernel(solver, NStates)
        probabilities = pauli.get_probabilities(solver, NStates)
        print("C++ probabilities:", probabilities)
        print("C++ kernel:\n", kernel)
    
    current = pauli.calculate_current(solver, 1)
    if verbosity > 0:
        print("C++ current:", current)
    
    pauli.cleanup(solver)
    return current
    
if __name__ == "__main__":
    # Define energy range
    bPrint = True
    nstep = 1
    eps = np.zeros((nstep,3))
    ts = np.linspace(0, 1, nstep)
    eps[:,0] = -10.0 + ts
    eps[:,1] = -10.01 + ts*2.0
    eps[:,2] = -10.02 + ts*0.3

    verbosity = 0
    
    # Run scan
    #qmeq_results = scan_QmeQ(eps, bPrint)
    #cpp_results  = scan_cpp(eps, bPrint)

    pauli = PauliSolver(verbosity=verbosity, bASAN=bASAN)

    qmeq_results = np.zeros(nstep)
    cpp_results = np.zeros(nstep)
    for i in range(nstep):
        epsi = eps[i]
        print(f"####### run python QmeQ Pauli {epsi}")
        qmeq_results[i] = run_QmeQ_solver(epsi[0], epsi[1], epsi[2])
        print(f"####### run C++ Pauli {epsi}")
        cpp_results[i] = run_cpp_solver(pauli, epsi[0], epsi[1], epsi[2])
        print(f"eps: {epsi} -> QmeQ: {qmeq_results[i]} C++: {cpp_results[i]}")
            
    plt.figure(figsize=(10, 6))
    plt.plot(ts, qmeq_results, 'o-b', label='QmeQ Pauli')
    plt.plot(ts, cpp_results,  'o:r', label='C++ Pauli')
    plt.xlabel('Onsite Energy (meV)')
    plt.ylabel('Current (nA)')
    plt.title('Solver Comparison for 1D Energy Scan')
    plt.legend()
    plt.grid(True)
    plt.show()