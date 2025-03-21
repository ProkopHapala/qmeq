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
from sys import path
path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')

import qmeq
from qmeq import indexing as qmqsi
#from qmeq.indexing import get_state_order
import traceback

# ==== Setup

# setup numpy print options to infinite line length
np.set_printoptions(linewidth=256, suppress=True)

import pauli_solver_lib as psl
from pauli_solver_lib import PauliSolver, calculate_state_energy, calculate_tunneling_amplitudes

# System parameters
NSingle = 3  # number of impurity states
NLeads  = 2   # number of leads

# Parameters (in meV) - with small perturbations to break degeneracy
eps1 = -10.0
eps2 = -10.01  # Slightly different
eps3 = -10.02  # Slightly different

t     = 0.0      # direct hopping
W     = 3.0     # inter-site coupling
VBias = 0.1  # bias voltage

# Lead parameters
muS    = 0.0    # substrate chemical potential
muT    = 0.0    # tip chemical potential
Temp   = 0.224 # temperature in meV
DBand  = 1000.0 # lead bandwidth
GammaS = 0.20 # coupling to substrate
GammaT = 0.05 # coupling to tip

# Tunneling amplitudes
VS = np.sqrt(GammaS/np.pi)  # substrate
VT = np.sqrt(GammaT/np.pi)  # tip

# Position-dependent coefficients
coeffE = 0.4
coeffT = 0.3


verbosity = 1

# ==== Functions

def build_hamiltonian(eps1, eps2, eps3, t, W):
    print("\n\n#### Building Hamiltonian: eps: ", [eps1, eps2, eps3], " t: ", t, " W: ", W)
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

def run_QmeQ_solver(Hsingle, Hcoulomb, mu_L, Temp_L, TLeads):  # Compare QmeQ and C++ solvers for same Hamiltonian
    # Run QmeQ solver
    if True:
        print( "\n\n" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "\n### Running QmeQ Pauli solver /home/prokop/git_SW/qmeq/qmeq/approach/base/pauli.py " )
        
        system = qmeq.Builder(NSingle, Hsingle, Hcoulomb, NLeads, TLeads, mu_L, Temp_L, DBand,   kerntype='Pauli', indexing='Lin', itype=0, symq=True,   solmethod='solve', mfreeq=0)
        system.appr.verbosity = verbosity  # Set verbosity after instance creation
        system.verbosity = verbosity
        system.solve()

        chargelst = system.si.chargelst
        state_order = qmqsi.get_state_order(chargelst); print("QmeQ state order:", state_order)
        state_occupancy = qmqsi.get_state_occupancy_strings(chargelst, NSingle); print("QmeQ state occupancy:", state_occupancy)
        
        res = {
            'current':        system.current[1],
            'energies':       system.Ea,
            'probabilities':  system.phi0,
            'kernel':         system.kern,
        }
        #print("\nQmeQ hsingle:", Hsingle)
        #print("QmeQ coulomb:", Hcoulomb)
        print("QmeQ energies:", system.Ea)
        print("QmeQ probabilities:", system.phi0)
        print("QmeQ kernel:\n", system.kern)

        return res

def run_cpp_solver(TLeads): 

    # Run C++ solver
    if True:
        print( "\n\n" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "\n### Running C++ solver /home/prokop/git_SW/qmeq/cpp/pauli_solver.hpp \n" )
        
        NStates = 2**NSingle
        
        lead_mu              = np.array([muS, muT + VBias])
        lead_temp            = np.array([Temp, Temp])
        lead_gamma           = np.array([GammaS, GammaT])
        pauli                = PauliSolver( verbosity=verbosity, bASAN=bASAN )

        TLeads_            = np.zeros( (NLeads, NStates) )
        for k,v in TLeads.items():
            TLeads_[k[0], k[1]] = v
        Hsingle_ = np.zeros( (NSingle, NSingle) )
        for k,v in Hsingle.items():
            Hsingle_[k[0], k[1]] = v

        print("\nHsingle:");         print(Hsingle_)
        print("\nTLeads:");         print(TLeads_)

        #state_order = [0, 1, 2, 4, 3, 5, 6, 7]
        state_order = [0, 4, 2, 6, 1, 5, 3, 7]
        state_order = np.array(state_order, dtype=np.int32)

        solver               = pauli.create_pauli_solver_new( NStates, NLeads, Hsingle_, W, TLeads_, lead_mu, lead_temp, lead_gamma, state_order, verbosity)

        energies = pauli.get_energies(solver, NStates)

        pauli.solve(solver)
        kernel               = pauli.get_kernel(solver, NStates)
        probabilities        = pauli.get_probabilities(solver, NStates)
        currents             = [pauli.calculate_current(solver, lead) for lead in range(NLeads)]
        pauli.cleanup(solver)
        
        res = {
            'current':       currents[1],
            'energies':      energies,
            'probabilities': probabilities,
            'kernel':        kernel,
        }

        #print("\nQmeQ hsingle:", Hsingle)
        #print("QmeQ coulomb:", Hcoulomb)
        print("C++ energies:", energies)
        print("C++ probabilities:", probabilities)
        print("C++ kernel:\n", kernel)

    
    return res

def test_cpp_kernel_with_numpy_solver(cpp_res):  # Verify C++ kernel by solving with NumPy
    """Test solving the C++ kernel matrix with NumPy's least squares solver"""
    print("\n\n#### Testing C++ kernel matrix with NumPy least squares solver:")
    
    # Get the kernel matrix from C++ result
    kernel = cpp_res['kernel'].copy()
    
    # Create the RHS vector [1, 0, 0, 0, 0, 0, 0, 0]
    n = kernel.shape[0]
    rhs = np.zeros(n)
    rhs[0] = 1.0
    
    # Replace first row with all ones (normalization condition)
    kernel[0, :] = 1.0
    
    # Print the kernel and RHS for verification
    print("\nC++ kernel matrix with normalization row:")
    np.set_printoptions(precision=15)
    print(kernel)
    print("\nRHS vector:")
    print(rhs)
    
    # Try least squares solver
    solution_lstsq, residuals, rank, s = np.linalg.lstsq(kernel, rhs, rcond=None)
    print("\nSolution from NumPy linalg.lstsq:")
    print(solution_lstsq)
    if len(residuals) > 0:
        print(f"Residuals: {residuals}")
    print(f"Rank: {rank}")
    print(f"Singular values: {s}")
    
    # Normalize the solution
    sum_prob = np.sum(solution_lstsq)
    if abs(sum_prob) > 1e-10:
        solution_normalized = solution_lstsq / sum_prob
    else:
        solution_normalized = solution_lstsq
    
    print("\nNormalized solution:")
    print(solution_normalized)
    
    # Compare with C++ probabilities
    print("\nC++ probabilities:")
    print(cpp_res['probabilities'])
    
    # Reset print options
    np.set_printoptions(precision=5, suppress=True)
    
    return solution_lstsq

def compare_results(qmeq_res, cpp_res, tol=1e-8, bPrintSame=True):  # Compare results from both solvers
    """Compare results from both solvers"""
    print( "\n\n" )
    print( "######################################################################" )
    print( "######################################################################" )
    print("#### Comparing QmeQ vs C++ results:")

    diff = np.max(np.abs(qmeq_res['energies'] - cpp_res['energies']))
    if diff > tol or bPrintSame:
        print("Energies diff :", diff)
        print("Energies QmeQ :", qmeq_res['energies'])
        print("Energies C++  :", cpp_res['energies'])
    else:
        print(f"Energies:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['probabilities'] - cpp_res['probabilities']))
    if diff > tol or bPrintSame:
        print("Probabilities diff :", diff)
        print("Probabilities QmeQ :", qmeq_res['probabilities'])
        print("Probabilities C++  :", cpp_res['probabilities'])
    else:
        print(f"Probabilities:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['kernel'] - cpp_res['kernel']))
    if diff > tol or bPrintSame:
        print("Kernel:   diff:", diff)
        print("Kernel QmeQ:\n", qmeq_res['kernel'])
        print("Kernel C++:\n", cpp_res['kernel'])
    else:
        print(f"Kernel:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['current'] - cpp_res['current']))
    if diff > tol or bPrintSame:
        print("Current diff :", diff)
        print("Current QmeQ :", qmeq_res['current'])
        print("Current C++  :", cpp_res['current'])
        print("Relative diff:", abs(qmeq_res['current'] - cpp_res['current'])/abs(qmeq_res['current']))
    else:
        print(f"Current:   OK (diff({diff}) < tol({tol}))")

# ==== Main

if __name__ == "__main__":

    mu_L, Temp_L, TLeads = build_leads(muS, muT, Temp, VS, VT, coeffT, VBias)
    Hsingle,  Hcoulomb = build_hamiltonian(eps1, eps2, eps3, t, W)
    print( "hsingle:\n", Hsingle)
    print( "coulomb:\n", Hcoulomb)

    cpp_res  = run_cpp_solver(TLeads)
    qmeq_res = run_QmeQ_solver(Hsingle, Hcoulomb, mu_L, Temp_L, TLeads)
    compare_results(qmeq_res, cpp_res)