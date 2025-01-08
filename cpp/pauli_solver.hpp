#pragma once

#include <vector>
#include <cstring>
#include <cstdio>
#include <cmath>
#include "gauss_solver.hpp"

// Constants should be defined in meV units
const double PI = 3.14159265358979323846;
const double HBAR = 0.6582119;  // Reduced Planck constant in meV*ps
const double KB = 0.08617333;   // Boltzmann constant in meV/K

// Lead parameters
struct LeadParams {
    double mu;    // Chemical potential
    double temp;  // Temperature
    double gamma; // Coupling strength
};

// Parameters for the solver
struct SolverParams {
    int nstates;  // Number of states
    int nleads;   // Number of leads
    std::vector<double> energies;  // State energies
    std::vector<LeadParams> leads; // Lead parameters
    std::vector<std::vector<std::vector<double>>> coupling; // Coupling matrix elements
};

class PauliSolver {
public:
    SolverParams params;
    double* kernel;         // Kernel matrix
    double* rhs;           // Right-hand side vector
    double* probabilities; // State probabilities
    double* pauli_factors; // Pauli factors for transitions
    int verbosity;        // Verbosity level for debugging
    std::vector<std::vector<int>> states_by_charge;  // States organized by charge number, like Python's statesdm

    // Count number of electrons in a state
    int count_electrons(int state) {
        return __builtin_popcount(state);
    }

    // Print matrix in numpy style
    void print_matrix(const double* mat, int rows, int cols, const char* label = nullptr) {
        if(label) printf("\n%s:\n", label);
        printf("[");
        for(int i = 0; i < rows; i++) {
            printf("[");
            for(int j = 0; j < cols; j++) {
                printf("%7.3f", mat[i * cols + j]);
                if (j < cols-1) printf(" ");
            }
            printf("]");
            if (i < rows-1) printf("\n ");
        }
        printf("]\n");
    }

    // Print vector in numpy style
    void print_vector(const double* vec, int size, const char* label = nullptr) {
        if(label) printf("\n%s:\n", label);
        printf("[");
        for(int i = 0; i < size; i++) {
            printf("%7.3f", vec[i]);
            if (i < size-1) printf(" ");
        }
        printf("]\n");
    }

    // Print Pauli factors array
    void print_pauli_factors(const char* label = nullptr) {
        if(label) printf("\n%s:\n", label);
        const int n = params.nstates;
        printf("Shape: [%d leads, %d states, %d states, 2 directions]\n", params.nleads, n, n);
        for(int l = 0; l < params.nleads; l++) {
            printf("\nLead %d:\n", l);
            for(int i = 0; i < n; i++) {
                printf(" [");
                for(int j = 0; j < n; j++) {
                    int idx = l * n * n * 2 + i * n * 2 + j * 2;
                    printf("[%7.3f %7.3f]", pauli_factors[idx], pauli_factors[idx + 1]);
                    if (j < n-1) printf(" ");
                }
                printf("]\n");
            }
        }
    }

    // Calculate Fermi function for given energy difference and lead parameters
    double fermi_func(double energy_diff, double mu, double temp) {
        return 1.0/(1.0 + exp((energy_diff - mu)/temp));
    }

    // Initialize states by charge number
    void init_states_by_charge() {
        const int n = params.nstates;
        int max_charge = 0;
        // Find maximum charge number
        for(int i = 0; i < n; i++) {
            max_charge = std::max(max_charge, count_electrons(i));
        }
        // Initialize the vector with empty vectors
        states_by_charge.resize(max_charge + 1);
        // Fill in states for each charge
        for(int i = 0; i < n; i++) {
            int charge = count_electrons(i);
            states_by_charge[charge].push_back(i);
        }
    }

    // Generate Pauli factors for transitions between states
    void generate_fct() {
        if(verbosity > 0) printf("\nDEBUG: generate_fct() Calculating Pauli factors...\n");
        
        const int n = params.nstates;
        memset(pauli_factors, 0, params.nleads * n * n * 2 * sizeof(double));
        
        // Make sure states are organized by charge
        if(states_by_charge.empty()) {
            init_states_by_charge();
        }
        
        // Iterate through charge states (like Python's implementation)
        for(int charge = 0; charge < states_by_charge.size() - 1; charge++) {
            int next_charge = charge + 1;
            
            // Iterate through states in current and next charge state
            for(int c : states_by_charge[next_charge]) {
                for(int b : states_by_charge[charge]) {
                    double energy_diff = params.energies[c] - params.energies[b];
                    
                    // For each lead
                    for(int l = 0; l < params.nleads; l++) {
                        const int idx = l * n * n * 2 + c * n * 2 + b * 2;
                        
                        // Calculate coupling strength
                        double coupling = params.coupling[l][b][c] * params.coupling[l][c][b];
                        
                        // Include lead parameters
                        const LeadParams& lead = params.leads[l];
                        double fermi = fermi_func(energy_diff, lead.mu, lead.temp);
                        
                        // Store factors for both directions
                        // Note: Python's func_pauli multiplies by 2π
                        pauli_factors[idx + 0] = coupling * lead.gamma * fermi * 2 * PI;         // Forward
                        pauli_factors[idx + 1] = coupling * lead.gamma * (1.0 - fermi) * 2 * PI; // Backward
                        
                        if(verbosity > 0) printf("DEBUG: generate_fct() l:%d i:%d j:%d E_diff:%.6f coupling:%.6f fermi:%.6f factors:[%.6f, %.6f]\n", 
                            l, c, b, energy_diff, coupling, fermi, pauli_factors[idx + 0], pauli_factors[idx + 1]);
                    }
                }
            }
        }
        
        if(verbosity > 0) {
            print_pauli_factors("Phase 1 - Pauli Factors");
        }
    }

    // Generate coupling terms for a specific state
    void generate_coupling_terms(int state) {
        if(verbosity > 0) printf("\nDEBUG: generate_coupling_terms() state:%d\n", state);
        
        const int n = params.nstates;
        int state_elec = count_electrons(state);
        double diagonal = 0.0;
        
        // Handle transitions from lower charge states (a → b)
        for(int other = 0; other < n; other++) {
            int other_elec = count_electrons(other);
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
            int other_elec = count_electrons(other);
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
        double* rhs = new double[n];
        std::fill(rhs, rhs + n - 1, 0.0);
        rhs[n - 1] = 1.0;
        
        // Solve the system using GaussSolver
        GaussSolver::solve(kern_copy, rhs, probabilities, n);
        
        delete[] kern_copy;
        delete[] rhs;
        
        if(verbosity > 0) {
            print_vector(probabilities, n, "Probabilities after solve");
        }
    }

    PauliSolver(const SolverParams& p, int verb = 0) : params(p), verbosity(verb) {
        const int n = params.nstates;
        kernel = new double[n * n];
        rhs = new double[n];
        probabilities = new double[n];
        pauli_factors = new double[params.nleads * n * n * 2];
    }

    ~PauliSolver() {
        delete[] kernel;
        delete[] rhs;
        delete[] probabilities;
        delete[] pauli_factors;
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
                int i_elec = count_electrons(i);
                int j_elec = count_electrons(j);
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
    const double* get_probabilities() const { return probabilities; }
    const double* get_rhs() const { return rhs; }
    const double* get_pauli_factors() const { return pauli_factors; }
};
