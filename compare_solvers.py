#!/usr/bin/env python3

import numpy as np
from sys import path
path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')
import qmeq

from pauli_solver_lib import PauliSolver, calculate_state_energy, calculate_tunneling_amplitudes

def run_qmeq_solver(params):
    """Run QmeQ solver with given parameters"""
    NSingle = params['NSingle']
    NLeads = params['NLeads']
    eps1, eps2, eps3 = params['eps1'], params['eps2'], params['eps3']
    t = params['t']
    W = params['W']
    VBias = params['VBias']
    muS, muT = params['muS'], params['muT']
    Temp = params['Temp']
    DBand = params['DBand']
    VS, VT = params['VS'], params['VT']
    coeffE, coeffT = params['coeffE'], params['coeffT']
    
    # One-particle Hamiltonian
    H1p = {(0,0): eps1-coeffE*VBias, (0,1): t, (0,2): t,
           (1,1): eps2, (1,2): t,
           (2,2): eps3}
    
    # Two-particle Hamiltonian: inter-site coupling
    H2p = {(0,1,1,0): W,
           (1,2,2,1): W,
           (0,2,2,0): W}
    
    # Leads: substrate (S) and scanning tip (T)
    mu_L = {0: muS, 1: muT + VBias}
    Temp_L = {0: Temp, 1: Temp}
    
    # Coupling between leads (1st number) and impurities (2nd number)
    TLeads = {(0,0): VS,         # S <-- 1
              (0,1): VS,         # S <-- 2
              (0,2): VS,         # S <-- 3
              (1,0): VT,         # T <-- 1
              (1,1): coeffT*VT,  # T <-- 2
              (1,2): coeffT*VT}  # T <-- 3
    
    # Create and solve system
    system = qmeq.Builder(NSingle, H1p, H2p, NLeads, TLeads, mu_L, Temp_L, DBand,
                         kerntype='Pauli', indexing='Lin', itype=0, symq=True,
                         solmethod='lsqr', mfreeq=0)
    system.solve()
    
    return {
        'current': system.current[1],
        'energies': system.Ea,
        'probabilities': system.phi0,
        'kernel': system.kern,
    }

def run_cpp_solver(params):
    """Run C++ solver with given parameters"""
    NSingle = params['NSingle']
    NStates = 2**NSingle
    NLeads = params['NLeads']
    eps1, eps2, eps3 = params['eps1'], params['eps2'], params['eps3']
    W = params['W']
    VS, VT = params['VS'], params['VT']
    muS, muT = params['muS'], params['muT']
    VBias = params['VBias']
    Temp = params['Temp']
    GammaS, GammaT = params['GammaS'], params['GammaT']
    coeffT = params['coeffT']
    
    # Calculate energies for each state
    energies = np.array([
        calculate_state_energy(i, NSingle, eps1, eps2, eps3, W) 
        for i in range(NStates)
    ])
    
    # Calculate tunneling amplitudes
    tunneling_amplitudes = calculate_tunneling_amplitudes(
        NLeads, NStates, NSingle, VS, VT, coeffT
    )
    
    # Lead parameters
    lead_mu = np.array([muS, muT + VBias])
    lead_temp = np.array([Temp, Temp])
    lead_gamma = np.array([GammaS, GammaT])
    
    # Create solver and run
    pauli = PauliSolver()
    solver = pauli.create_solver(
        NStates, NLeads, energies, tunneling_amplitudes,
        lead_mu, lead_temp, lead_gamma
    )
    
    # Solve and get results
    pauli.solve(solver)
    kernel = pauli.get_kernel(solver, NStates)
    probabilities = pauli.get_probabilities(solver, NStates)
    currents = [pauli.calculate_current(solver, lead) for lead in range(NLeads)]
    
    # Cleanup
    pauli.cleanup(solver)
    
    return {
        'current': currents[1],  # tip current
        'energies': energies,
        'probabilities': probabilities,
        'kernel': kernel,
    }

def compare_results(qmeq_res, cpp_res):
    """Compare results from both solvers"""
    print("Comparing QmeQ vs C++ results:")
    print("\nEnergies:")
    print("QmeQ:", qmeq_res['energies'])
    print("C++: ", cpp_res['energies'])
    print("Max diff:", np.max(np.abs(qmeq_res['energies'] - cpp_res['energies'])))
    
    print("\nProbabilities:")
    print("QmeQ:", qmeq_res['probabilities'])
    print("C++: ", cpp_res['probabilities'])
    print("Max diff:", np.max(np.abs(qmeq_res['probabilities'] - cpp_res['probabilities'])))
    
    print("\nKernel:")
    print("QmeQ:", qmeq_res['kernel'])
    print("C++: ", cpp_res['kernel'])
    print("Max diff:", np.max(np.abs(qmeq_res['kernel'] - cpp_res['kernel'])))
    
    print("\nCurrent:")
    print("QmeQ:", qmeq_res['current'])
    print("C++: ", cpp_res['current'])
    print("Relative diff:", abs(qmeq_res['current'] - cpp_res['current'])/abs(qmeq_res['current']))

def main():
    # System parameters
    params = {
        'NSingle': 3,  # number of impurity states
        'NLeads': 2,   # number of leads
        
        # Parameters (in meV)
        'eps1': -10.0,
        'eps2': -10.0,
        'eps3': -10.0,
        't': 0.0,      # direct hopping
        'W': 20.0,     # inter-site coupling
        'VBias': 0.0,  # bias voltage
        
        # Lead parameters
        'muS': 0.0,    # substrate chemical potential
        'muT': 0.0,    # tip chemical potential
        'Temp': 0.224, # temperature in meV
        'DBand': 1000.0, # lead bandwidth
        'GammaS': 0.20, # coupling to substrate
        'GammaT': 0.05, # coupling to tip
        
        # Tunneling amplitudes
        'VS': np.sqrt(0.20/np.pi),  # substrate
        'VT': np.sqrt(0.05/np.pi),  # tip
        
        # Position-dependent coefficients
        'coeffE': 0.4,
        'coeffT': 0.3,
    }
    
    # Run both solvers
    qmeq_res = run_qmeq_solver(params)
    cpp_res = run_cpp_solver(params)
    
    # Compare results
    compare_results(qmeq_res, cpp_res)

if __name__ == "__main__":
    main()
