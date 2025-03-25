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

import pauli_lib as psl
from pauli_lib import PauliSolver

# setup numpy print options to infinite line length
np.set_printoptions(linewidth=256, suppress=True)

# Constants
NSingle = 3  # number of impurity states
NLeads = 2   # number of leads

# Parameters (in meV) - with small perturbations to break degeneracy
eps1 = -10.0
eps2 = -10.01  # Slightly different
eps3 = -10.02  # Slightly different

t = 0.0      # direct hopping
W = 3.0      # inter-site coupling (matching compare_solvers.py)
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
    if verbosity > 0:
        print( "\n\n" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "\n### Running QmeQ Pauli solver /home/prokop/git_SW/qmeq/qmeq/approach/base/pauli.py " )
    
    mu_L, Temp_L, TLeads = build_leads(muS, muT, Temp, VS, VT, coeffT, VBias)
    Hsingle, Hcoulomb    = build_hamiltonian(eps1, eps2, eps3, t, W)
    
    try:
        config.verbosity = verbosity
        system = qmeq.Builder(NSingle, Hsingle, Hcoulomb, NLeads, TLeads, mu_L, Temp_L, DBand,  kerntype='Pauli', indexing='Lin', itype=0, symq=True, solmethod='solve', mfreeq=0)
        system.appr.verbosity = verbosity  # Set verbosity after instance creation
        system.verbosity = verbosity
        if verbosity > 0:
            print( "type(system).__name__ ", type(system).__name__)
        system.solve()
    except Exception as e:
        print(f"Error running QmeQ solver: {e}")
        return None
    
    chargelst = system.si.chargelst
    state_order = qmqsi.get_state_order(chargelst); 
    state_occupancy = qmqsi.get_state_occupancy_strings(chargelst, NSingle); 
    
    if verbosity > 0:
        print("QmeQ state order:", state_order)
        print("QmeQ state occupancy:", state_occupancy)
        print("QmeQ energies:", system.Ea)
        print("QmeQ probabilities:", system.phi0)
        print("QmeQ kernel:\n", system.kern)
        print("QmeQ current:", system.current[1])
    
    # Create a result dictionary for comparison
    qmeq_res = {
        'current': system.current[1],
        'energies': system.Ea,
        'probabilities': system.phi0,
        'kernel': system.kern,
        'pauli_factors': system.appr.paulifct,
        'leads': {
            'mu': system.leads.mulst,
            'temp': system.leads.tlst,
            'gamma': system.leads.dlst[:,0],
            'Tba': system.leads.Tba.real
        }
    }
    
    return qmeq_res

def prepare_cpp_inputs(eps1, eps2, eps3):
    """Prepare inputs for C++ solver"""
    
    mu_L, Temp_L, TLeads = build_leads(muS, muT, Temp, VS, VT, coeffT, VBias)
    Hsingle, Hcoulomb    = build_hamiltonian(eps1, eps2, eps3, t, W)
    lead_mu    = np.array([muS    , muT + VBias ])
    lead_temp  = np.array([Temp   , Temp        ])
    lead_gamma = np.array([GammaS , GammaT      ])
    
    # Convert TLeads dictionary to matrix form
    TLeads_ = np.zeros((NLeads, NSingle))
    for k, v in TLeads.items():
        TLeads_[k[0], k[1]] = v
    
    # Convert Hsingle dictionary to matrix form
    Hsingle_ = np.zeros((NSingle, NSingle))
    for k, v in Hsingle.items():
        Hsingle_[k[0], k[1]] = v
    if verbosity > 0:
        print("\nHsingle:");        print(Hsingle_)
        print("\nTLeads:");         print(TLeads_)
    return Hsingle_, TLeads_, lead_mu, lead_temp, lead_gamma

def prepare_cpp_inputs_efficient(eps1, eps2, eps3):
    """Efficiently prepare inputs for C++ solver using numpy arrays directly"""
    # Single particle Hamiltonian
    Hsingle = np.array([
        [eps1, t, 0],
        [t, eps2, t], 
        [0, t, eps3]])
    
    # Leads
    lead_mu    = np.array([muS,    muT + VBias ] )
    lead_temp  = np.array([Temp,   Temp        ] )
    lead_gamma = np.array([GammaS, GammaT      ] )  
    # Lead Tunneling matrix
    TLeads = np.array([
        [VS, VS        , VS       ],
        [VT, VT*coeffT , VT*coeffT]])
    
    return Hsingle, TLeads, lead_mu, lead_temp, lead_gamma

def prepare_leads_cpp():
    """Prepare static inputs that don't change with eps"""
    # Leads
    lead_mu = np.array([muS, muT + VBias])
    lead_temp = np.array([Temp, Temp])
    lead_gamma = np.array([GammaS, GammaT])
    # Lead Tunneling matrix
    TLeads = np.array([
        [VS, VS, VS],
        [VT, VT*coeffT, VT*coeffT]
    ])    
    return TLeads, lead_mu, lead_temp, lead_gamma

def prepare_hsinglecpp(eps1, eps2, eps3):
    """Prepare dynamic inputs that change with eps"""
    # Single particle Hamiltonian
    Hsingle = np.array([
        [eps1, t, 0],
        [t, eps2, t],
        [0, t, eps3]
    ])
    return Hsingle

def run_cpp_solver(pauli, eps1, eps2, eps3):
    """Run C++ solver with the given parameters"""
    if verbosity > 0:
        print( "\n\n" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "\n### Running C++ solver /home/prokop/git_SW/qmeq/cpp/pauli_solver.hpp \n" )
    
    #Hsingle_, TLeads_, lead_mu, lead_temp, lead_gamma = prepare_cpp_inputs(eps1, eps2, eps3)
    #Hsingle_, TLeads_, lead_mu, lead_temp, lead_gamma = prepare_cpp_inputs_efficient(eps1, eps2, eps3)
    
    # --- prepare static inputs - this does not change when we change eps
    state_order = [0, 4, 2, 6, 1, 5, 3, 7]
    state_order = np.array(state_order, dtype=np.int32)
    TLeads_, lead_mu, lead_temp, lead_gamma = prepare_leads_cpp()
    NStates = 2**NSingle

    # --- prepare dynamic inputs - this changes when we change eps
    Hsingle_ = prepare_hsinglecpp(eps1, eps2, eps3)
    
    # Create and run solver
    solver = pauli.create_pauli_solver_new(NStates, NLeads, Hsingle_, W, TLeads_, lead_mu, lead_temp, lead_gamma, state_order, verbosity)
    
    # Get energies before solving
    energies = pauli.get_energies(solver, NStates)
    if verbosity > 0:
        print("C++ energies:", energies)
    
    pauli.solve(solver)
    
    # Get detailed results for comparison
    kernel         = pauli.get_kernel(solver, NStates)
    probabilities = pauli.get_probabilities(solver, NStates)
    currents      = [pauli.calculate_current(solver, lead) for lead in range(NLeads)]
    Tba           = pauli.get_coupling(solver, NLeads, NStates)
    pauli_factors = pauli.get_pauli_factors(solver, NLeads, NStates)
    
    if verbosity > 0:
        print("C++ probabilities:", probabilities)
        print("C++ kernel:\n", kernel)
        print("C++ current:", currents[1])
    
    # Create a result dictionary for comparison
    cpp_res = {
        'current': currents[1],
        'energies': energies,
        'probabilities': probabilities,
        'kernel': kernel,
        'pauli_factors': pauli_factors,
        'leads': {
            'mu': lead_mu,
            'temp': lead_temp,
            'gamma': lead_gamma,
            'Tba': Tba
        }
    }
    
    pauli.cleanup(solver)
    return cpp_res
    

# Define a comparison function like in compare_solvers.py
def compare_results(qmeq_res, cpp_res, tol=1e-8, bPrintSame=True):
    print("\n\n### Comparing QmeQ and C++ results")
    
    # Compare current
    qmeq_current = qmeq_res['current']
    cpp_current  = cpp_res['current']
    diff_current = abs(qmeq_current - cpp_current)
    
    # Compare energies
    qmeq_energies = qmeq_res['energies']
    cpp_energies = cpp_res['energies']
    
    print("\nEnergies:")
    for i, (qe, ce) in enumerate(zip(qmeq_energies, cpp_energies)):
        diff = abs(qe - ce)
        if diff > tol or bPrintSame:
            print(f"  State {i}: QmeQ={qe}, C++={ce}, Diff={diff}")
    
    # Compare probabilities
    qmeq_probs = qmeq_res['probabilities']
    cpp_probs  = cpp_res['probabilities']
    
    print("\nProbabilities:")
    for i, (qp, cp) in enumerate(zip(qmeq_probs, cpp_probs)):
        diff = abs(qp - cp)
        if diff > tol or bPrintSame:
            print(f"  State {i}: QmeQ={qp}, C++={cp}, Diff={diff}")

if __name__ == "__main__":
    print( "\n\n" )
    print( "##################################################################################" )
    print( "##################################################################################" )
    print( "### compare_scan_1D.py Compare QmeQ vs C++ Pauli solvers for 1D array of energies" )
    print( "##################################################################################" )
    print( "##################################################################################" )
    

    # Use exact parameters from compare_solvers.py
    bPrint = True
    nstep = 5
    # Instead of generating energy range, use exact values from compare_solvers.py
    eps = np.zeros((nstep,3))
    ts = np.linspace(0, 1.0, nstep)
    eps[:,0] = eps1 + ts   # eps1 from compare_solvers.py
    eps[:,1] = eps2 + ts  # eps2 from compare_solvers.py
    eps[:,2] = eps3 + ts  # eps3 from compare_solvers.py

    verbosity = 0  # Match compare_solvers.py verbosity
    
    pauli = PauliSolver(verbosity=verbosity, bASAN=bASAN)

    Iqmeq = np.zeros(nstep)
    Icpp  = np.zeros(nstep)
    

    for i in range(nstep):
        epsi = eps[i]
        if verbosity > 0: print(f"\n####### compare_scan_1D.py loop [{i}] epsi: {epsi}")
        qmeq_res = run_QmeQ_solver(epsi[0], epsi[1], epsi[2])
        cpp_res  = run_cpp_solver(pauli, epsi[0], epsi[1], epsi[2])
        print(f"Eps[{i}] {epsi} Current: QmeQ: {qmeq_res['current']} C++: {cpp_res['current']} | Diff: {abs(qmeq_res['current'] - cpp_res['current'])}")
        #compare_results(qmeq_res, cpp_res, tol=1e-8, bPrintSame=True)
        Iqmeq[i] = qmeq_res['current']
        Icpp[i]  = cpp_res['current']

    plt.figure(figsize=(10, 6))
    plt.plot(ts, Iqmeq, 'o-b', label='QmeQ Pauli')
    plt.plot(ts, Icpp,  'o:r', label='C++ Pauli')
    plt.xlabel('Onsite Energy (meV)')
    plt.ylabel('Current (nA)')
    plt.title('Solver Comparison for 1D Energy Scan')
    plt.legend()
    plt.grid(True)
    plt.show()