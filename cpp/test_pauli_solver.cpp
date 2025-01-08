#include "pauli_solver.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <memory>
#include <fstream>
#include <bitset>

// Parameters matching spinless_trimer1.py
const int NSingle = 3;
const int NStates = 1 << NSingle;  // 2^NSingle
const int NLeads = 2;

// System parameters (in meV)
const double eps1 = -10.0;
const double eps2 = -10.0;
const double eps3 = -10.0;
const double t = 0.0;    // direct hopping
const double W = 20.0;   // inter-site coupling
const double VBias = 0.0;  // We'll test with zero bias first

// Lead parameters
const double muS = 0.0;    // substrate chemical potential
const double muT = 0.0;    // tip chemical potential
const double Temp = 0.224; // temperature in meV
const double GammaS = 0.20;  // coupling to substrate
const double GammaT = 0.05;  // coupling to tip
const double VS = std::sqrt(GammaS/M_PI);  // tunneling amplitude substrate
const double VT = std::sqrt(GammaT/M_PI);  // tunneling amplitude tip

// Helper function to count electrons in a state
int count_electrons(int state) {
    return __builtin_popcount(state);
}

// Helper function to calculate state energy
double calculate_state_energy(int state) {
    double energy = 0.0;
    // Count number of electrons and their positions
    for(int i = 0; i < NSingle; i++) {
        if(state & (1 << i)) {
            energy += (i == 0) ? eps1 : ((i == 1) ? eps2 : eps3);
        }
    }
    // Add inter-site coupling W for adjacent occupied sites
    for(int i = 0; i < NSingle-1; i++) {
        if((state & (1 << i)) && (state & (1 << (i+1)))) {
            energy += W;
        }
    }
    return energy;
}

// Helper function to check if states differ by one electron at a specific site
bool is_valid_transition(int state1, int state2, int site) {
    int diff = state1 ^ state2;
    return (diff == (1 << site));
}

int main() {
    // Create system parameters
    SolverParams params;
    params.nstates = NStates;
    params.nleads = NLeads;
    
    // Allocate and fill energies
    params.energies.resize(NStates);
    for(int i = 0; i < NStates; i++) {
        params.energies[i] = calculate_state_energy(i);
    }
    
    // Allocate and fill coupling matrix
    params.coupling.resize(NLeads);
    for(int l = 0; l < NLeads; l++) {
        params.coupling[l].resize(NStates);
        for(int i = 0; i < NStates; i++) {
            params.coupling[l][i].resize(NStates);
        }
    }
    
    // For each lead and each pair of states
    for(int lead = 0; lead < NLeads; lead++) {
        double v = (lead == 0) ? VS : VT;  // Choose coupling based on lead
        
        for(int i = 0; i < NStates; i++) {
            for(int j = 0; j < NStates; j++) {
                // States must differ by exactly one electron
                int diff = count_electrons(i) - count_electrons(j);
                if(std::abs(diff) != 1) continue;
                
                // Check if transition is valid (only one site changes)
                bool valid = false;
                for(int site = 0; site < NSingle; site++) {
                    if(is_valid_transition(i, j, site)) {
                        valid = true;
                        double v_mod = v;
                        // Apply position-dependent coupling for tip
                        if(lead == 1) {  // Tip
                            double coeff = (site == 0) ? 1.0 : 0.3;  // coeffT = 0.3 from Python
                            v_mod *= coeff;
                        }
                        params.coupling[lead][i][j] = v_mod;
                        break;
                    }
                }
            }
        }
    }
    
    // Initialize lead parameters
    params.leads.resize(NLeads);
    
    // Substrate lead
    params.leads[0].mu = muS;
    params.leads[0].temp = Temp;
    params.leads[0].gamma = GammaS * 2.0 * M_PI;  // Convert from QmeQ's Γ/(2π) to Γ
    
    // Tip lead
    params.leads[1].mu = muT;
    params.leads[1].temp = Temp;
    params.leads[1].gamma = GammaT * 2.0 * M_PI;  // Convert from QmeQ's Γ/(2π) to Γ
    
    // Create solver
    PauliSolver solver(params, 1);  // verbosity = 1
    
    // Solve
    solver.solve();
    
    return 0;
}
