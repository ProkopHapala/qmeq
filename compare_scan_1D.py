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
    # Initialize QmeQ solver
    hsingle = {(0,0): 0.0, (0,1): t, (0,2): t,
               (1,1): 0.0, (1,2): t,
               (2,2): 0.0}
    coulomb = {(0,1,1,0): W, (1,2,2,1): W, (0,2,2,0): W}
    mu_L = {0: muS, 1: muT + VBias}
    Temp_L = {0: Temp, 1: Temp}

    qmeq_system = qmeq.Builder(NSingle, hsingle, coulomb, NLeads, TLeads, mu_L, Temp_L, DBand, kerntype='Pauli', indexing='Lin', itype=0, symq=True, solmethod='solve', mfreeq=0)

    return qmeq_system

def initialize_cpp_solver():
    """Initialize C++ solver once"""
    pauli = PauliSolver()
    NStates = 2**NSingle
    
    # Create constant parts of the solver
    tunneling_amplitudes = calculate_tunneling_amplitudes(NLeads, NStates, NSingle, TLeads)
    lead_mu = np.array([muS, muT + VBias])
    lead_temp = np.array([Temp, Temp])
    lead_gamma = np.array([GammaS, GammaT])
    
    # Create solver instance with dummy energies
    solver = pauli.create_solver(NStates, NLeads, np.zeros(NStates), tunneling_amplitudes, lead_mu, lead_temp, lead_gamma)
    
    return pauli, solver, NStates

def run_qmeq_solver(eps1, eps2, eps3, qmeq_system):
    """Run QmeQ solver with updated Hamiltonian"""
    # Update QmeQ Hamiltonian
    qmeq_system.hsingle[(0,0)] = eps1-coeffE*VBias
    qmeq_system.hsingle[(1,1)] = eps2
    qmeq_system.hsingle[(2,2)] = eps3
    qmeq_system.solve()
    return qmeq_system.current[1]

def run_cpp_solver(eps1, eps2, eps3, pauli, solver, NStates):
    """Run C++ solver with updated energies"""
    # Calculate new energies
    energies = np.array([calculate_state_energy(i, NSingle, eps1, eps2, eps3, W, VBias=VBias, coeffE=coeffE, t=t) for i in range(NStates)])
    
    # Recreate solver with new energies
    tunneling_amplitudes = calculate_tunneling_amplitudes(NLeads, NStates, NSingle, TLeads)
    lead_mu    = np.array([muS, muT + VBias])
    lead_temp  = np.array([Temp, Temp])
    lead_gamma = np.array([GammaS, GammaT])
    
    # Create new solver instance
    solver = pauli.create_solver(NStates, NLeads, energies, tunneling_amplitudes, lead_mu, lead_temp, lead_gamma)
    
    # Solve and get current
    pauli.solve(solver)
    return pauli.calculate_current(solver, 1)

def scan_energies(eps_range, solver='both'):
    """Scan through a range of onsite energies"""
    results = {'cpp': [], 'qmeq': [], 'eps': eps_range}
    
    # Initialize solvers once
    if solver in ['qmeq', 'both']:
        qmeq_system = initialize_qmeq_solver()
    if solver in ['cpp', 'both']:
        pauli, cpp_solver, NStates = initialize_cpp_solver()
    
    for eps in eps_range:
        if solver in ['qmeq', 'both']:
            qmeq_current = run_qmeq_solver(eps, eps, eps, qmeq_system)
            results['qmeq'].append(qmeq_current)
        
        if solver in ['cpp', 'both']:
            cpp_current = run_cpp_solver(eps, eps, eps, pauli, cpp_solver, NStates)
            results['cpp'].append(cpp_current)
    
    return results

def plot_results(results):
    """Plot comparison of solver results"""
    plt.figure(figsize=(10, 6))
    
    if results['qmeq']:
        plt.plot(results['eps'], results['qmeq'], 'b-', label='QmeQ Pauli')
    if results['cpp']:
        plt.plot(results['eps'], results['cpp'], 'r--', label='C++ Pauli')
    
    plt.xlabel('Onsite Energy (meV)')
    plt.ylabel('Current (nA)')
    plt.title('Solver Comparison for 1D Energy Scan')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    # Define energy range
    eps_range = np.linspace(-20, 20, 50)
    
    # Run scan
    results = scan_energies(eps_range, solver='both')
    
    # Plot results
    plot_results(results)
