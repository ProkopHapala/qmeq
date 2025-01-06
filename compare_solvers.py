#!/usr/bin/env python3

import numpy as np
from sys import path
path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')
import qmeq

# setup numpy print options to infinite line length
np.set_printoptions(linewidth=256, suppress=True)
#np.set_printoptions(linewidth=256, suppress=False)

from pauli_solver_lib import PauliSolver, calculate_state_energy, calculate_tunneling_amplitudes

def run_solvers( bRunQmeq=True, bRunCpp=True ):

    # System parameters
    NSingle = 3  # number of impurity states
    NLeads  = 2   # number of leads

    # Parameters (in meV)
    eps1 = eps2 = eps3 = -10.0
    t     = 0.0      # direct hopping
    W     = 20.0     # inter-site coupling
    VBias = 0.0  # bias voltage

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
    
    # One-particle Hamiltonian
    hsingle = {(0,0): eps1-coeffE*VBias, (0,1): t, (0,2): t,
               (1,1): eps2, (1,2): t,
               (2,2): eps3}
    
    # Two-particle Hamiltonian: inter-site coupling
    coulomb = {(0,1,1,0): W,
               (1,2,2,1): W,
               (0,2,2,0): W}
    
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


    verbosity = 1
    
    # Run QmeQ solver
    if bRunQmeq:
        print( "\n\n### Running QmeQ Pauli solver /home/prokop/git_SW/qmeq/qmeq/approach/base/pauli.py " )
        system = qmeq.Builder(NSingle, hsingle, coulomb, NLeads, TLeads, mu_L, Temp_L, DBand,   kerntype='Pauli', indexing='Lin', itype=0, symq=True,   solmethod='lsqr', mfreeq=0)
        
        system.verbosity = verbosity
        system.solve()
        qmeq_res = {
            'current': system.current[1],
            'energies': system.Ea,
            'probabilities': system.phi0,
            'kernel': system.kern,
        }
        print("\nQmeQ hsingle:", hsingle)
        print("QmeQ coulomb:", coulomb)
        print("QmeQ states:", [bin(i)[2:].zfill(NSingle) for i in range(2**NSingle)])
        print("QmeQ energies:", system.Ea)

    # Run C++ solver
    if bRunCpp:
        print( "\n\n### Running C++ solver /home/prokop/git_SW/qmeq/cpp/pauli_solver.hpp " )
        NStates = 2**NSingle
        energies = np.array([calculate_state_energy(i, NSingle, eps1, eps2, eps3, W, VBias=VBias, coeffE=coeffE, t=t) 
                            for i in range(NStates)])
        tunneling_amplitudes = calculate_tunneling_amplitudes(  NLeads, NStates, NSingle, VS, VT, coeffT )
        lead_mu              = np.array([muS, muT + VBias])
        lead_temp            = np.array([Temp, Temp])
        lead_gamma           = np.array([GammaS, GammaT])
        pauli                = PauliSolver()
        solver               = pauli.create_solver(   NStates, NLeads, energies, tunneling_amplitudes,  lead_mu, lead_temp, lead_gamma, verbosity)
        pauli.solve(solver)
        kernel               = pauli.get_kernel(solver, NStates)
        probabilities        = pauli.get_probabilities(solver, NStates)
        currents             = [pauli.calculate_current(solver, lead) for lead in range(NLeads)]
        pauli.cleanup(solver)
        
        cpp_res = {
            'current': currents[1],
            'energies': energies,
            'probabilities': probabilities,
            'kernel': kernel,
        }
        print("\nC++ states:", [bin(i)[2:].zfill(NSingle) for i in range(NStates)])
        print("C++ energies:", energies)
    
    return qmeq_res, cpp_res

def compare_results(qmeq_res, cpp_res, tol=1e-8):
    """Compare results from both solvers"""
    print("\n\n#### Comparing QmeQ vs C++ results:")

    diff = np.max(np.abs(qmeq_res['energies'] - cpp_res['energies']))
    if diff > tol:
        print("Energies:   diff:", diff)
        print("QmeQ:", qmeq_res['energies'])
        print("C++: ", cpp_res['energies'])
    else:
        print(f"Energies:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['probabilities'] - cpp_res['probabilities']))
    if diff > tol:
        print("Probabilities:   diff:", diff)
        print("QmeQ:", qmeq_res['probabilities'])
        print("C++: ", cpp_res['probabilities'])
    else:
        print(f"Probabilities:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['kernel'] - cpp_res['kernel']))
    if diff > tol:
        print("Kernel:   diff:", diff)
        print("QmeQ:", qmeq_res['kernel'])
        print("C++: ", cpp_res['kernel'])
    else:
        print(f"Kernel:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['current'] - cpp_res['current']))
    if diff > tol:
        print("Current:   diff:", diff)
        print("QmeQ:", qmeq_res['current'])
        print("C++: ", cpp_res['current'])
        print("Relative diff:", abs(qmeq_res['current'] - cpp_res['current'])/abs(qmeq_res['current']))
    else:
        print(f"Current:   OK (diff({diff}) < tol({tol}))")

def main():
    qmeq_res, cpp_res = run_solvers()
    compare_results(qmeq_res, cpp_res)

if __name__ == "__main__":
    main()
