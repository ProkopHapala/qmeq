#!/usr/bin/env python3

import numpy as np
from sys import path
path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')

import qmeq
import traceback

# ==== Setup

# setup numpy print options to infinite line length
np.set_printoptions(linewidth=256, suppress=True)

from pauli_solver_lib import PauliSolver, calculate_state_energy, calculate_tunneling_amplitudes

# System parameters
NSingle = 3  # number of impurity states
NLeads  = 2   # number of leads

# Parameters (in meV) - with small perturbations to break degeneracy
eps1 = -10.0
eps2 = -10.01  # Slightly different
eps3 = -10.02  # Slightly different

t     = 0.0      # direct hopping
W     = 20.0     # inter-site coupling
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

# ==== Functions

def build_hamiltonian(eps1, eps2, eps3, t, W):
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

def run_QmeQ_solver(Hsingle, Hcoulomb, mu_L, Temp_L, TLeads, verbosity=4):  # Compare QmeQ and C++ solvers for same Hamiltonian
    # Run QmeQ solver
    if True:
        print( "\n\n" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "\n### Running QmeQ Pauli solver /home/prokop/git_SW/qmeq/qmeq/approach/base/pauli.py " )
        
        try:
            system = qmeq.Builder(NSingle, Hsingle, Hcoulomb, NLeads, TLeads, mu_L, Temp_L, DBand,   kerntype='Pauli', indexing='Lin', itype=0, symq=True,   solmethod='solve', mfreeq=0)
            system.appr.verbosity = verbosity  # Set verbosity after instance creation
            system.verbosity = verbosity
            system.solve()

        except Exception as e:
            print(f"Error in QmeQ solver: {e}")
            traceback.print_exc()

        res = {
            'current': system.current[1],
            'energies': system.Ea,
            'probabilities': system.phi0,
            'kernel': system.kern,
        }
        #print("\nQmeQ hsingle:", Hsingle)
        #print("QmeQ coulomb:", Hcoulomb)
        print("QmeQ states:", [bin(i)[2:].zfill(NSingle) for i in range(2**NSingle)])
        print("QmeQ energies:", system.Ea)
        print("QmeQ probabilities:", system.phi0)
        print("QmeQ kernel:", system.kern)

        return res

def run_cpp_solver(Hsingle, Hcoulomb, mu_L, Temp_L, TLeads, verbosity=4): 

    # Run C++ solver
    if True:
        print( "\n\n" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "\n### Running C++ solver /home/prokop/git_SW/qmeq/cpp/pauli_solver.hpp \n" )
        # print("\n System parameters:")
        # print(f"NSingle={NSingle}, NLeads={NLeads}")
        # print(f"eps1={eps1}, eps2={eps2}, eps3={eps3}")
        # print(f"t={t}, W={W}, VBias={VBias}")
        # print(f"GammaS={GammaS}, GammaT={GammaT}")
        # print(f"VS={VS}, VT={VT}")
        # print(f"coeffE={coeffE}, coeffT={coeffT}")
        
        NStates = 2**NSingle
        energies = np.array([calculate_state_energy(i, NSingle, eps1, eps2, eps3, W, VBias=VBias, coeffE=coeffE, t=t) for i in range(NStates)])
        tunneling_amplitudes = calculate_tunneling_amplitudes(NLeads, NStates, NSingle, TLeads)
        
        lead_mu              = np.array([muS, muT + VBias])
        lead_temp            = np.array([Temp, Temp])
        lead_gamma           = np.array([GammaS, GammaT])
        pauli                = PauliSolver( verbosity=verbosity )
        solver               = pauli.create_solver(   NStates, NLeads, energies, tunneling_amplitudes,  lead_mu, lead_temp, lead_gamma, verbosity)
        pauli.solve(solver)
        kernel               = pauli.get_kernel(solver, NStates)
        probabilities        = pauli.get_probabilities(solver, NStates)
        currents             = [pauli.calculate_current(solver, lead) for lead in range(NLeads)]
        pauli.cleanup(solver)
        
        res = {
            'current': currents[1],
            'energies': energies,
            'probabilities': probabilities,
            'kernel': kernel,
        }

        #print("\nQmeQ hsingle:", Hsingle)
        #print("QmeQ coulomb:", Hcoulomb)
        print("\nC++ states:", [bin(i)[2:].zfill(NSingle) for i in range(NStates)])
        print("C++ energies:", energies)
        print("C++ probabilities:", probabilities)
        print("C++ kernel:", kernel)

    
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

def compare_results(qmeq_res, cpp_res, tol=1e-8):  # Compare results from both solvers
    """Compare results from both solvers"""
    print( "\n\n" )
    print( "######################################################################" )
    print( "######################################################################" )
    print("#### Comparing QmeQ vs C++ results:")

    diff = np.max(np.abs(qmeq_res['energies'] - cpp_res['energies']))
    if diff > tol:
        print("Energies diff :", diff)
        print("Energies QmeQ :", qmeq_res['energies'])
        print("Energies C++  :", cpp_res['energies'])
    else:
        print(f"Energies:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['probabilities'] - cpp_res['probabilities']))
    if diff > tol:
        print("Probabilities diff :", diff)
        print("Probabilities QmeQ :", qmeq_res['probabilities'])
        print("Probabilities C++  :", cpp_res['probabilities'])
    else:
        print(f"Probabilities:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['kernel'] - cpp_res['kernel']))
    if diff > tol:
        print("Kernel:   diff:", diff)
        print("Kernel QmeQ:\n", qmeq_res['kernel'])
        print("Kernel C++:\n", cpp_res['kernel'])
    else:
        print(f"Kernel:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['current'] - cpp_res['current']))
    if diff > tol:
        print("Current diff :", diff)
        print("Current QmeQ :", qmeq_res['current'])
        print("Current C++  :", cpp_res['current'])
        print("Relative diff:", abs(qmeq_res['current'] - cpp_res['current'])/abs(qmeq_res['current']))
    else:
        print(f"Current:   OK (diff({diff}) < tol({tol}))")

# ==== Main

if __name__ == "__main__":

    mu_L, Temp_L, TLeads = build_leads(muS, muT, Temp, VS, VT, coeffT, VBias)

    # -- Symmetric Hamiltonian
    #print("\n\n#### Testing with symmetric Hamiltonian (with degeneracies):")
    # Hsingle, Hcoulomb = build_hamiltonian(eps1, eps1, eps1, t, W)
    # qmeq_res = run_QmeQ_solver(Hsingle, Hcoulomb, mu_L, Temp_L, TLeads)
    # cpp_res  = run_cpp_solver(Hsingle, Hcoulomb, mu_L, Temp_L, TLeads)
    # compare_results(qmeq_res, cpp_res)
    
    # -- Asymmetric Hamiltonian
    #print("\n\n#### Testing with perturbed Hamiltonian to break degeneracy:")
    #eps1, eps2, eps3 = eps1, eps2, eps3
    Hsingle,  Hcoulomb = build_hamiltonian(eps1, eps2, eps3, t, W)
    cpp_res  = run_cpp_solver (Hsingle, Hcoulomb, mu_L, Temp_L, TLeads)
    qmeq_res = run_QmeQ_solver(Hsingle, Hcoulomb, mu_L, Temp_L, TLeads)
    compare_results(qmeq_res, cpp_res)