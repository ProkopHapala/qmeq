#ifndef PAULI_SOLVER_HPP
#define PAULI_SOLVER_HPP

#include <cmath>
#include <cstring>
#include <cassert>

// Constants
const double PI = 3.14159265358979323846;
const double HBAR = 1.0545718e-34;  // Reduced Planck constant
const double KB = 1.380649e-23;     // Boltzmann constant

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

class PauliSolver {
private:
    SystemParams params;
    double* kernel;     // Kernel matrix L [nstates x nstates]
    double* rhs;        // Right-hand side vector b [nstates]
    double* probabilities; // State probabilities p [nstates]

    // Calculate transition rate W_ji from state j to i
    double calculate_rate(int i, int j, int lead_idx) {
        const int n = params.nstates;
        double energy_diff = params.energies[i] - params.energies[j];
        double tunneling = params.tunneling_amplitudes[lead_idx * n * n + j * n + i];
        double coupling = tunneling * tunneling;  // |⟨i|H_tunneling|j⟩|^2
        
        // Include lead parameters
        const LeadParams& lead = params.leads[lead_idx];
        double fermi = fermi_func(energy_diff, lead.mu, lead.temp);
        
        // Rate according to Fermi's golden rule
        return (2.0 * PI / HBAR) * coupling * lead.gamma * fermi;
    }

    // Build the kernel matrix L and right-hand side vector b
    void build_system() {
        const int n = params.nstates;
        memset(kernel, 0, n * n * sizeof(double));
        memset(rhs, 0, n * sizeof(double));

        // Build kernel matrix L
        for(int i = 0; i < n; i++) {
            double diagonal_sum = 0.0;
            
            // Off-diagonal elements (transition rates)
            for(int j = 0; j < n; j++) {
                if(i == j) continue;
                
                double total_rate = 0.0;
                for(int l = 0; l < params.nleads; l++) {
                    total_rate += calculate_rate(i, j, l);
                }
                
                kernel[i * n + j] = total_rate;  // W_ji
                diagonal_sum += total_rate;
            }
            
            // Diagonal elements (negative sum of outgoing rates)
            kernel[i * n + i] = -diagonal_sum;
        }

        // Replace last row with normalization condition
        for(int j = 0; j < n; j++) {
            kernel[(n-1) * n + j] = 1.0;
        }
        rhs[n-1] = 1.0;  // Σp_i = 1
    }

    // Solve the system Lp = b using Gaussian elimination with partial pivoting
    void solve_system() {
        const int n = params.nstates;
        double* aug_matrix = new double[n * (n + 1)];  // Augmented matrix [L|b]
        
        // Create augmented matrix
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                aug_matrix[i * (n + 1) + j] = kernel[i * n + j];
            }
            aug_matrix[i * (n + 1) + n] = rhs[i];
        }

        // Gaussian elimination with partial pivoting
        for(int k = 0; k < n - 1; k++) {
            // Find pivot
            int pivot_row = k;
            double pivot_val = fabs(aug_matrix[k * (n + 1) + k]);
            for(int i = k + 1; i < n; i++) {
                double val = fabs(aug_matrix[i * (n + 1) + k]);
                if(val > pivot_val) {
                    pivot_val = val;
                    pivot_row = i;
                }
            }

            // Swap rows if necessary
            if(pivot_row != k) {
                for(int j = k; j <= n; j++) {
                    double temp = aug_matrix[k * (n + 1) + j];
                    aug_matrix[k * (n + 1) + j] = aug_matrix[pivot_row * (n + 1) + j];
                    aug_matrix[pivot_row * (n + 1) + j] = temp;
                }
            }

            // Eliminate column
            for(int i = k + 1; i < n; i++) {
                double factor = aug_matrix[i * (n + 1) + k] / aug_matrix[k * (n + 1) + k];
                for(int j = k; j <= n; j++) {
                    aug_matrix[i * (n + 1) + j] -= factor * aug_matrix[k * (n + 1) + j];
                }
            }
        }

        // Back substitution
        for(int i = n - 1; i >= 0; i--) {
            probabilities[i] = aug_matrix[i * (n + 1) + n];
            for(int j = i + 1; j < n; j++) {
                probabilities[i] -= aug_matrix[i * (n + 1) + j] * probabilities[j];
            }
            probabilities[i] /= aug_matrix[i * (n + 1) + i];
        }

        delete[] aug_matrix;
    }

public:
    PauliSolver(const SystemParams& p) : params(p) {
        const int n = params.nstates;
        kernel = new double[n * n];
        rhs = new double[n];
        probabilities = new double[n];
    }

    ~PauliSolver() {
        delete[] kernel;
        delete[] rhs;
        delete[] probabilities;
    }

    // Solve the Pauli master equation and return state probabilities
    const double* solve() {
        build_system();
        solve_system();
        return probabilities;
    }

    // Calculate current through a specific lead
    double calculate_current(int lead_idx) {
        double current = 0.0;
        const int n = params.nstates;
        
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(i == j) continue;
                // Current = e * Σ_ij (W_ji^(l) p_j - W_ij^(l) p_i)
                double rate_ji = calculate_rate(i, j, lead_idx);
                double rate_ij = calculate_rate(j, i, lead_idx);
                current += rate_ji * probabilities[j] - rate_ij * probabilities[i];
            }
        }
        
        return current;  // in units of electron charge
    }
};

#endif // PAULI_SOLVER_HPP
