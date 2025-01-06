#ifndef PAULI_SOLVER_HPP
#define PAULI_SOLVER_HPP

#include <cmath>
#include <cstring>
#include <cassert>
#include <cstdio>
#include "gauss_solver.hpp"

// Constants should be defined in meV units
const double PI = 3.14159265358979323846;
const double HBAR = 0.6582119;  // Reduced Planck constant in meV*ps
const double KB = 0.08617333;   // Boltzmann constant in meV/K

// Structure to hold lead parameters
struct LeadParams {
    double mu;      // Chemical potential
    double temp;    // Temperature
    double gamma;   // Coupling strength (wide-band limit)
};

// Structure to hold system parameters
struct SystemParams {
    int nstates;    // Total number of many-body states
    int nleads;     // Number of leads
    double* energies;  // Array of state energies [nstates]
    double* tunneling_amplitudes;  // Array of tunneling amplitudes [nleads][nstates][nstates]
    LeadParams* leads;  // Array of lead parameters [nleads]
};

// Fermi-Dirac distribution function
inline double fermi_func(double energy, double mu, double temp) {
    if(temp == 0.0) return (energy <= mu) ? 1.0 : 0.0;
    return 1.0 / (1.0 + exp((energy - mu) / temp));
}

void print_matrix(double* matrix, int nrows, int ncols) {
    for(int i = 0; i < nrows; i++) {
        for(int j = 0; j < ncols; j++) {
            printf("%.6f ", matrix[i * ncols + j]);
        }
        printf("\n");
    }
}

class PauliSolver { public:
    SystemParams params;
    double* kernel;        // Kernel matrix L [nstates x nstates]
    double* rhs;          // Right-hand side vector b [nstates]
    double* probabilities; // State probabilities p [nstates]
    double* pauli_factors; // Transition factors [nleads][nstates][nstates][2]
    int verbosity = 0;
    
    // Calculate Pauli factors for a transition
    void generate_fct() {
        if(verbosity > 0) printf("\nDEBUG: generate_fct() Calculating Pauli factors...\n");
        const int n = params.nstates;
        const int nl = params.nleads;
        
        // For each lead and pair of states
        for(int l = 0; l < nl; l++) {
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    int idx = l * n * n * 2 + i * n * 2 + j * 2;
                    double energy_diff = params.energies[i] - params.energies[j];
                    double tunneling = params.tunneling_amplitudes[l * n * n + j * n + i];
                    double coupling = tunneling * tunneling;  // |⟨i|H_tunneling|j⟩|^2
                    
                    // Include lead parameters
                    const LeadParams& lead = params.leads[l];
                    double fermi = fermi_func(energy_diff, lead.mu, lead.temp);
                    
                    // Store factors for both directions (like QmeQ's paulifct)
                    pauli_factors[idx + 0] = coupling * lead.gamma * fermi * (2.0 * PI) / HBAR;        // Forward
                    pauli_factors[idx + 1] = coupling * lead.gamma * (1.0 - fermi) * (2.0 * PI) / HBAR; // Backward
                    
                    if(verbosity > 0) printf("DEBUG: generate_fct() l:%d i:%d j:%d E_diff:%.6f coupling:%.6f fermi:%.6f factors:[%.6f, %.6f]\n",  l, i, j, energy_diff, coupling, fermi, pauli_factors[idx + 0], pauli_factors[idx + 1]);
                }
            }
        }
    }

    // Set kernel matrix elements for a specific state
    void generate_coupling_terms(int state) {
        if(verbosity > 0) printf("\nDEBUG: generate_coupling_terms() state:%d\n", state);
        const int n = params.nstates;
        const int nl = params.nleads;
        
        double diagonal_sum = 0.0;
        
        // Handle transitions to other states
        for(int other = 0; other < n; other++) {
            if(other == state) continue;
            
            int state_elec = count_electrons(state);
            int other_elec = count_electrons(other);
            if(abs(state_elec - other_elec) != 1) continue; // Only single electron transitions
            
            double total_rate = 0.0;
            for(int l = 0; l < nl; l++) {
                int idx = l * n * n * 2 + state * n * 2 + other * 2;
                // Use appropriate factor based on electron number difference
                total_rate += (state_elec > other_elec) ? 
                             pauli_factors[idx + 1] :  // Electron leaving
                             pauli_factors[idx + 0];   // Electron entering
            }
            
            kernel[state * n + other] = total_rate;
            diagonal_sum += total_rate;
            
            if(verbosity > 0) printf("DEBUG: generate_coupling_terms() state:%d other:%d rate:%.6f\n",   state, other, total_rate);
        }
        
        // Set diagonal element
        kernel[state * n + state] = -diagonal_sum;
        if(verbosity > 0) printf("DEBUG: generate_coupling_terms() state:%d diagonal:%.6f\n",   state, -diagonal_sum);
    }

    

    // Build the full kernel matrix
    void generate_kern() {
        if(verbosity > 0) printf("\nDEBUG: generate_kern() Building kernel matrix...\n");
        const int n = params.nstates;
        
        // Initialize arrays
        memset(kernel, 0, n * n * sizeof(double));
        memset(rhs, 0, n * sizeof(double));
        
        // Generate kernel matrix elements for each state
        for(int state = 0; state < n-1; state++) {
            generate_coupling_terms(state);
        }
        
        // Replace last row with normalization condition
        for(int j = 0; j < n; j++) {
            kernel[(n-1) * n + j] = 1.0;
        }
        rhs[n-1] = 1.0;  // Σp_i = 1
        
        if(verbosity > 0){ 
            printf("DEBUG: generate_kern() Kernel matrix:\n");
            print_matrix( kernel, n, n);
        }
    }

    // Count number of electrons in a state
    int count_electrons(int state) {
        int count = 0;
        while(state) {
            count += state & 1;
            state >>= 1;
        }
        return count;
    }

public:
    PauliSolver(const SystemParams& p, int verbosity_=0 ) : params(p), verbosity(verbosity_) {
        const int n = params.nstates;
        const int nl = params.nleads;
        kernel = new double[n * n];
        rhs = new double[n];
        probabilities = new double[n];
        pauli_factors = new double[nl * n * n * 2]; // [nleads][nstates][nstates][2]
    }

    ~PauliSolver() {
        delete[] kernel;
        delete[] rhs;
        delete[] probabilities;
        delete[] pauli_factors;
    }

    // Solve the Pauli master equation
    void solve() {
        generate_fct();
        generate_kern();
        
        const int n = params.nstates;
        GaussSolver::solve(kernel, rhs, probabilities, n);
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
                if(verbosity > 0) printf("DEBUG: generate_current() i:%d j:%d rate:%.6f prob:%.6f contrib:%.6f\n",  i, j, rate, probabilities[j], rate * probabilities[j] * (i_elec > j_elec ? -1.0 : 1.0));
            }
        }
        return current;
    }

    // Getter methods
    const double* get_kernel() const { return kernel; }
    const double* get_probabilities() const { return probabilities; }
    int get_nstates() const { return params.nstates; }
};

#endif // PAULI_SOLVER_HPP
