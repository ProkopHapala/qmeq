#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <bitset>
#include "print_utils.hpp"
#include "gauss_solver.hpp"

// Constants should be defined in meV units
const double PI = 3.14159265358979323846;
const double kB = 8.617333262e-2;  // Boltzmann constant in meV/K

struct LeadParams {
    double mu;    // Chemical potential
    double temp;  // Temperature
    double gamma; // Coupling strength
};

// System parameters
struct SystemParams {
    // Basic system parameters
    int nsingle;  // Number of single-particle states
    int nstates;  // Number of many-body states (2^nsingle)
    int nleads;   // Number of leads
    
    // Single-particle parameters
    double eps1;  // Site energies
    double eps2;
    double eps3;
    
    double t;     // Hopping parameter
    double W;     // Inter-site coupling
    double VBias; // Bias voltage
    
    double GammaS; // Coupling to substrate
    double GammaT; // Coupling to tip
    double coeffE; // Energy coefficient
    double coeffT; // Tunneling coefficient
    
    double muS;   // Substrate chemical potential
    double muT;   // Tip chemical potential
    double Temp;  // Temperature
    
    std::vector<LeadParams> leads; // Lead parameters
    
    // Calculated arrays
    std::vector<double> energies;
    std::vector<std::vector<std::vector<double>>> coupling;
    
    SystemParams() {} // Default constructor
    
    SystemParams(int nsingle_, double eps1_, double eps2_, double eps3_,
                double t_, double W_, double VBias_,
                double GammaS_, double GammaT_, double coeffE_, double coeffT_,
                double muS_, double muT_, double Temp_)
        : nsingle(nsingle_), nstates(1 << nsingle_), nleads(2),
          eps1(eps1_), eps2(eps2_), eps3(eps3_),
          t(t_), W(W_), VBias(VBias_),
          GammaS(GammaS_), GammaT(GammaT_), coeffE(coeffE_), coeffT(coeffT_),
          muS(muS_), muT(muT_), Temp(Temp_) {
        
        // Calculate VS and VT from GammaS/T and W
        //double VS = sqrt(GammaS * W / (2 * PI));
        //double VT = sqrt(GammaT * W / (2 * PI));
        
        leads.resize(nleads);
        leads[0] = {muS_, Temp_, GammaS_ * 2.0 * M_PI};
        leads[1] = {muT_, Temp_, GammaT_ * 2.0 * M_PI};
        
        // Initialize arrays
        energies.resize(nstates);
        coupling.resize(nleads);
        for(int l = 0; l < nleads; l++) {
            coupling[l].resize(nstates);
            for(int i = 0; i < nstates; i++) {
                coupling[l][i].resize(nstates, 0.0);
            }
        }
    }
};

class PauliSolver {
private:
    SystemParams params;
    int verbosity;
    std::vector<std::vector<int>> states_by_charge;
    std::vector<double> pauli_factors;
    std::vector<double> probabilities;
    double* kernel;
    double* rhs;
    
    // Get number of electrons in a state
    int get_state_charge(int state) {
        return __builtin_popcount(state);
    }

    // Get position of the different bit between two states
    int get_transition_site(int state1, int state2) {
        int diff = state1 ^ state2;
        if(__builtin_popcount(diff) != 1) return -1;
        return __builtin_ctz(diff);  // Count trailing zeros
    }

    // Calculate Fermi function for given energy difference and lead parameters
    double fermi_func(double energy_diff, double mu, double temp) {
        double x = (energy_diff - mu)/temp;
        if (x > 700.0) return 0.0;  // exp(x) would overflow
        if (x < -700.0) return 1.0; // exp(x) would underflow
        return 1.0/(1.0 + exp(x));
    }
    
    // Calculate energy of a given state
    double calculate_state_energy(int state) {
        double energy = 0.0;
        
        // Add site energies for occupied sites
        if(state & 1) energy += params.eps1;
        if(state & 2) energy += params.eps2;
        if(state & 4) energy += params.eps3;
        
        // Add interaction energy between electrons
        if(params.W != 0.0) {
            // Between sites 1-2
            if((state & 1) && (state & 2)) energy += params.W;
            // Between sites 2-3
            if((state & 2) && (state & 4)) energy += params.W;
            // Between sites 1-3
            if((state & 1) && (state & 4)) energy += params.W;
        }
        
        // Add hopping terms
        if(params.t != 0.0) {
            // Between sites 1-2
            if((state & 1) && !(state & 2)) energy += params.t;
            if(!(state & 1) && (state & 2)) energy += params.t;
            
            // Between sites 2-3
            if((state & 2) && !(state & 4)) energy += params.t * params.coeffT;
            if(!(state & 2) && (state & 4)) energy += params.t * params.coeffT;
        }
        
        // Add bias voltage contribution
        if(params.VBias != 0.0 && params.coeffE != 0.0) {
            double bias_factor = 0.0;
            if(state & 1) bias_factor += 0.0;        // Site 1
            if(state & 2) bias_factor += 0.5;        // Site 2
            if(state & 4) bias_factor += 1.0;        // Site 3
            energy += params.VBias * params.coeffE * bias_factor;
        }
        
        return energy;
    }
    
    // Calculate all state energies
    void calculate_energies() {
        if(verbosity > 0) printf("\nDEBUG: Calculating state energies...\n");
        
        for(int i = 0; i < params.nstates; i++) {
            params.energies[i] = calculate_state_energy(i);
            if(verbosity > 0) {
                printf("State %s: %.6f\n", 
                       std::bitset<3>(i).to_string().c_str(), 
                       params.energies[i]);
            }
        }
    }
    
    // Calculate tunneling amplitudes
    void calculate_tunneling_amplitudes() {
        if(verbosity > 0) {
            printf("\nDEBUG: Calculating tunneling amplitudes...\n");
            printf("coeffT=%.6f\n", params.coeffT);
        }
        
        // For each lead and pair of states
        for(int l = 0; l < params.nleads; l++) {
            double v = (l == 0) ? sqrt(params.leads[l].gamma / (2 * PI)) : sqrt(params.leads[l].gamma / (2 * PI));
            
            for(int i = 0; i < params.nstates; i++) {
                for(int j = 0; j < params.nstates; j++) {
                    // Count number of different bits (must be 1 for valid transition)
                    int diff = i ^ j;
                    if(__builtin_popcount(diff) != 1) continue;
                    
                    // Get position of different bit (site index)
                    int site = __builtin_ctz(diff);
                    
                    // Calculate amplitude based on lead and site
                    double coeff = (l == 0 || site == 0) ? 1.0 : params.coeffT;
                    double amplitude = v * coeff;
                    
                    if(verbosity > 0) {
                        printf("l:%d i:%d j:%d site:%d v:%.6f coeff:%.6f amplitude:%.6f\n",
                               l, i, j, site, v, coeff, amplitude);
                    }
                    
                    params.coupling[l][i][j] = amplitude;
                }
            }
        }
    }

    // Initialize states by charge number
    void init_states_by_charge() {
        if(verbosity > 0) printf("\nDEBUG: Grouping states by charge...\n");
        
        int max_charge = 0;
        for(int i = 0; i < params.nstates; i++) {
            max_charge = std::max(max_charge, get_state_charge(i));
        }
        
        states_by_charge.resize(max_charge + 1);
        for(int i = 0; i < params.nstates; i++) {
            int charge = get_state_charge(i);
            states_by_charge[charge].push_back(i);
            
            if(verbosity > 0) {
                printf("State %s: charge %d\n", 
                       std::bitset<3>(i).to_string().c_str(), 
                       charge);
            }
        }
    }

public:
    PauliSolver(const SystemParams& p, int verb = 0) : params(p), verbosity(verb) {
        if(verbosity > 0) {
            printf("\nDEBUG: System parameters:\n");
            printf("NSingle=%d, NLeads=%d\n", params.nsingle, params.nleads);
            printf("eps1=%.1f, eps2=%.1f, eps3=%.1f\n", params.eps1, params.eps2, params.eps3);
            printf("t=%.1f, W=%.1f, VBias=%.1f\n", params.t, params.W, params.VBias);
            printf("GammaS=%.3f, GammaT=%.3f\n", params.GammaS, params.GammaT);
            printf("coeffE=%.1f, coeffT=%.1f\n", params.coeffE, params.coeffT);
        }
        
        // Allocate arrays
        kernel = new double[params.nstates * params.nstates];
        rhs = new double[params.nstates];
        probabilities.resize(params.nstates);
        pauli_factors.resize(params.nleads * params.nstates * params.nstates * 2, 0.0);
        
        // Initialize everything
        calculate_energies();
        calculate_tunneling_amplitudes();
        init_states_by_charge();
    }

    ~PauliSolver() {
        delete[] kernel;
        delete[] rhs;
    }

    // Generate Pauli factors for transitions between states
    void generate_fct() {
        if(verbosity > 0) printf("\nDEBUG: generate_fct() Calculating Pauli factors...\n");
        
        const int n = params.nstates;
        pauli_factors.resize(params.nleads * n * n * 2, 0.0);
        
        // Make sure states are organized by charge
        if(states_by_charge.empty()) {
            init_states_by_charge();
        }
        
        // Iterate through charge states (like Python's implementation)
        for(int charge = 0; charge < (int) states_by_charge.size() - 1; charge++) {
            int next_charge = charge + 1;
            
            // Iterate through states in current and next charge state
            for(int c : states_by_charge[next_charge]) {
                for(int b : states_by_charge[charge]) {
                    double energy_diff = params.energies[c] - params.energies[b];
                    
                    // For each lead
                    for(int l = 0; l < params.nleads; l++) {
                        const int idx = l * n * n * 2 + c * n * 2 + b * 2;
                        
                        // Get the site that changed in this transition
                        int changed_site = get_transition_site(c, b);
                        
                        // Calculate coupling strength
                        double coupling = params.coupling[l][b][c] * params.coupling[l][c][b];
                        
                        // Apply position-dependent coupling for tip (lead 1)
                        if (l == 1 && changed_site > 0) {  // For tip and sites 1,2
                            coupling *= 0.09;  // (coeffT * coeffT) = 0.3 * 0.3
                        }
                        
                        // Include lead parameters
                        double mu = params.leads[l].mu;
                        double fermi = fermi_func(energy_diff, mu, params.leads[l].temp);
                        
                        // Store factors for both directions
                        // Note: Python's func_pauli multiplies by 2π and coupling already includes gamma/π
                        pauli_factors[idx + 0] = coupling * fermi * 2 * PI;         // Forward
                        pauli_factors[idx + 1] = coupling * (1.0 - fermi) * 2 * PI; // Backward
                        
                        if(verbosity > 0) printf("DEBUG: generate_fct() l:%d i:%d j:%d E_diff:%.6f coupling:%.6f fermi:%.6f factors:[%.6f, %.6f]\n", 
                            l, c, b, energy_diff, coupling, fermi, pauli_factors[idx + 0], pauli_factors[idx + 1]);
                    }
                }
            }
        }
        
        if(verbosity > 0) {
            printf("\nDEBUG: Pauli factors:\n");
            for(int l = 0; l < params.nleads; l++) {
                printf("Lead %d:\n", l);
                for(int dir = 0; dir < 2; dir++) {
                    printf("Direction %d:\n", dir);
                    print_matrix(&pauli_factors[l * params.nstates * params.nstates * 2 + dir], 
                               params.nstates, params.nstates);
                }
            }
        }
    }

    // Generate coupling terms for a specific state
    void generate_coupling_terms(int state) {
        if(verbosity > 0) printf("\nDEBUG: generate_coupling_terms() state:%d\n", state);
        
        const int n = params.nstates;
        int state_elec = get_state_charge(state);
        double diagonal = 0.0;
        
        // Handle transitions from lower charge states (a → b)
        for(int other = 0; other < n; other++) {
            int other_elec = get_state_charge(other);
            if(other_elec != state_elec - 1) continue;  // Only from lower charge states
            
            double rate = 0.0;
            for(int l = 0; l < params.nleads; l++) {
                int idx = l * n * n * 2 + other * n * 2 + state * 2;
                rate += pauli_factors[idx + 0];  // Electron entering
            }
            kernel[other * n + state] = rate;  // Off-diagonal term
            diagonal -= rate;  // Add to diagonal term
        }
        
        // Handle transitions to higher charge states (b → c)
        for(int other = 0; other < n; other++) {
            int other_elec = get_state_charge(other);
            if(other_elec != state_elec + 1) continue;  // Only to higher charge states
            
            double rate = 0.0;
            for(int l = 0; l < params.nleads; l++) {
                int idx = l * n * n * 2 + other * n * 2 + state * 2;
                rate += pauli_factors[idx + 1];  // Electron leaving
            }
            kernel[other * n + state] = rate;  // Off-diagonal term
            diagonal -= rate;  // Add to diagonal term
        }
        
        // Set diagonal term
        kernel[state * n + state] = diagonal;
        
        if(verbosity > 0) {
            printf("DEBUG: generate_coupling_terms() state:%d diagonal:%.6f\n", 
                state, diagonal);
            print_matrix(kernel, n, n, "Phase 2 - After processing state");
        }
    }

    // Generate kernel matrix
    void generate_kern() {
        if(verbosity > 0) printf("\nDEBUG: generate_kern() Building kernel matrix...\n");
        
        const int n = params.nstates;
        
        // Calculate Pauli factors first
        generate_fct();
        
        // Initialize kernel matrix to zero
        std::fill(kernel, kernel + n * n, 0.0);
        
        // Generate coupling terms for each state
        for(int state = 0; state < n; state++) {
            generate_coupling_terms(state);
        }
        
        // Replace last row with normalization condition
        for(int j = 0; j < n; j++) {
            kernel[n * (n-1) + j] = 1.0;  // Sum of probabilities = 1
        }
        
        if(verbosity > 0) {
            print_matrix(kernel, n, n, "Phase 2 - After normalization");
        }
    }

    // Solve the kernel matrix equation
    void solve_kern() {
        const int n = params.nstates;
        
        // Create a copy of kernel matrix since solve() modifies it
        double* kern_copy = new double[n * n];
        std::copy(kernel, kernel + n * n, kern_copy);
        
        // Set up RHS vector [0, 0, ..., 1]
        std::fill(rhs, rhs + n - 1, 0.0);
        rhs[n - 1] = 1.0;
        
        // Solve the system using GaussSolver
        GaussSolver::solve(kern_copy, rhs, &probabilities[0], n);
        
        delete[] kern_copy;
        
        if(verbosity > 0) {
            print_vector(&probabilities[0], n, "Probabilities after solve");
        }
    }

    // Solve the master equation
    void solve() {
        generate_fct();
        solve_kern();
    }

    // Calculate current through a specific lead
    double generate_current(int lead_idx) {
        if(verbosity > 0) printf("\nDEBUG: generate_current() lead:%d\n", lead_idx);
        
        const int n = params.nstates;
        double current = 0.0;
        
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                int i_elec = get_state_charge(i);
                int j_elec = get_state_charge(j);
                if(abs(i_elec - j_elec) != 1) continue;
                
                int idx = lead_idx * n * n * 2 + j * n * 2 + i * 2;
                double rate = (i_elec > j_elec) ? 
                             pauli_factors[idx + 1] :  // Electron leaving
                             pauli_factors[idx + 0];   // Electron entering
                
                current += rate * probabilities[j] * (i_elec > j_elec ? -1.0 : 1.0);
                if(verbosity > 0) printf("DEBUG: generate_current() i:%d j:%d rate:%.6f prob:%.6f contrib:%.6f\n",
                    i, j, rate, probabilities[j], rate * probabilities[j] * (i_elec > j_elec ? -1.0 : 1.0));
            }
        }
        return current;
    }

    // Getter methods
    const double* get_kernel() const { return kernel; }
    const double* get_probabilities() const { return &probabilities[0]; }
    const double* get_rhs() const { return rhs; }
    const double* get_pauli_factors() const { return &pauli_factors[0]; }
};
