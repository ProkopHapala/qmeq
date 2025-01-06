#!/usr/bin/env python3

import numpy as np
from pauli_solver_lib import PauliSolver, calculate_state_energy, calculate_tunneling_amplitudes

def main():
    # System parameters
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
    
    print("\nCurrents:")
    print(f"Substrate: {currents[0]}")
    print(f"Tip: {currents[1]}")
    
    # Cleanup
    pauli.cleanup(solver)

if __name__ == "__main__":
    main()
