#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from sys import path

import qmeq
from qmeq import config
from qmeq import indexing as qmqsi
from qmeq.config import verb_print_

import pauli_solver_lib as psl
from pauli_solver_lib import PauliSolver

# Constants
NSingle = 3  # number of impurity states
NLeads = 2   # number of leads

# Parameters (in meV)
t = 0.0      # direct hopping
W = 3.0      # inter-site coupling
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

# Scan parameters
eps1_range = np.linspace(-20, 20, 100)
eps2_range = np.linspace(-20, 20, 100)
eps3 = -10.02  # Fixed value

bQmeQ=True
bCpp=True
verbosity = 0

def build_leads(muS, muT, Temp, VS, VT, coeffT, VBias):
    """Build leads"""
    mu_L   = {0: muS, 1: muT + VBias}
    Temp_L = {0: Temp, 1: Temp}
    TLeads = {(0,0): VS,         # S <-- 1
              (0,1): VS,         # S <-- 2
              (0,2): VS,         # S <-- 3
              (1,0): VT,         # T <-- 1
              (1,1): coeffT*VT,  # T <-- 2
              (1,2): coeffT*VT}  # T <-- 3
    return mu_L, Temp_L, TLeads

def build_hamiltonian(eps1, eps2, eps3, t, W):
    """Build Hamiltonian"""
    Hsingle  = np.array([[eps1, t, 0], [t, eps2, t], [0, t, eps3]])
    Hcoulomb = np.array([[W, 0, 0], [0, W, 0], [0, 0, W]])
    return Hsingle, Hcoulomb

def scan_qmeq(eps1_range, eps2_range, eps3, t, W):
    """Perform 2D scan using QmeQ solver"""
    currents = np.zeros((len(eps1_range), len(eps2_range)))
    mu_L, Temp_L, TLeads = build_leads(muS, muT, Temp, VS, VT, coeffT, VBias)
    Hsingle, Hcoulomb = build_hamiltonian(eps1_range[0], eps2_range[0], eps3, t, W)
    system = qmeq.Builder(NSingle, Hsingle, Hcoulomb, NLeads, TLeads, mu_L, Temp_L, DBand, kerntype='Pauli', indexing='Lin', itype=0, symq=True, solmethod='solve', mfreeq=0)
    # Perform scan
    for i, eps1 in enumerate(eps1_range):
        for j, eps2 in enumerate(eps2_range):
            Hsingle, Hcoulomb = build_hamiltonian(eps1, eps2, eps3, t, W)
            system.update_hamiltonian(Hsingle, Hcoulomb)
            system.solve()
            currents[i,j] = system.current[1]
    return currents


def scan_cpp(eps1_range, eps2_range, eps3, t, W):
    """Perform 2D scan using C++ solver"""
    currents = np.zeros((len(eps1_range), len(eps2_range)))
    pauli = PauliSolver()
    mu_L, Temp_L, TLeads = build_leads(muS, muT, Temp, VS, VT, coeffT, VBias)
    state_order = [0, 4, 2, 6, 1, 5, 3, 7]  # Example state order
    for i, eps1 in enumerate(eps1_range):
        for j, eps2 in enumerate(eps2_range):
            Hsingle, Hcoulomb = build_hamiltonian(eps1, eps2, eps3, t, W)
            solver = pauli.create_pauli_solver_new(2**NSingle, NLeads, Hsingle, W, TLeads, mu_L, Temp_L, [GammaS, GammaT], state_order, verbosity)
            pauli.solve(solver)
            currents[i,j] = pauli.calculate_current(solver, 1)
    
    return currents


def plot_results(eps1_range, eps2_range, currents, label='Current', cmap='viridis', title=None):
    """Plot results from both solvers"""
    plt.figure(figsize=(12,5))
    plt.subplot(111)
    plt.imshow(currents, extent=[eps2_range[0], eps2_range[-1], eps1_range[0], eps1_range[-1]], origin='lower', aspect='auto', cmap=cmap)
    plt.colorbar(label=label)
    plt.xlabel('eps2 (meV)')
    plt.ylabel('eps1 (meV)')
    if title:
        plt.title(title)

if __name__ == "__main__":
    
    # if bQmeQ:
    #     qmeq_currents = scan_qmeq(eps1_range, eps2_range, eps3, t, W)
    #     plot_results(eps1_range, eps2_range, qmeq_currents, label='Current (QmeQ)')

    if bCpp:
        cpp_currents = scan_cpp(eps1_range, eps2_range, eps3, t, W)
        plot_results(eps1_range, eps2_range, cpp_currents, label='Current (C++)')

    plt.tight_layout()
    plt.show()
