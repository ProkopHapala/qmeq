#ifndef PAULI_SOLVER_HPP
#define PAULI_SOLVER_HPP

#include <cmath>
#include <cstring>
#include <cassert>
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

class PauliSolver {
private:
    SystemParams params;
    double* kernel;        // Kernel matrix L [nstates x nstates]
    double* rhs;           // Right-hand side vector b [nstates]
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
        // Note: QmeQ uses Γ = 2π|t|²ρ, where ρ=1/2π in wide-band limit
        // So our coupling already includes the 2π factor
        // Note: QmeQ multiplies by 2π in func_pauli
        return coupling * lead.gamma * fermi * (2.0 * PI) / HBAR;
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
        GaussSolver::solve(kernel, rhs, probabilities, n);
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

    // Solve the Pauli master equation
    void solve() {
        build_system();
        solve_system();
    }

    // Getter methods
    const double* get_kernel() const { return kernel; }
    const double* get_probabilities() const { return probabilities; }
    int get_nstates() const { return params.nstates; }

    // Calculate current through a specific lead
    double calculate_current(int lead_idx) {
        const int n = params.nstates;
        double current = 0.0;
        
        // For each pair of states
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                // Get tunneling amplitude and rate
                double tunneling = params.tunneling_amplitudes[lead_idx * n * n + j * n + i];
                if(tunneling == 0.0) continue;  // Skip if no tunneling
                
                // Calculate current contribution
                // Note: QmeQ uses I = e ∑_ij W_ij p_j where W_ij is the rate from j to i
                // The sign is determined by whether i has more electrons than j
                int i_elec = count_electrons(i);
                int j_elec = count_electrons(j);
                if(i_elec == j_elec) continue;  // Skip if same number of electrons
                
                double rate = calculate_rate(i, j, lead_idx);
                current += rate * probabilities[j] * (i_elec > j_elec ? 1.0 : -1.0);
            }
        }
        return current;
    }
};

#endif // PAULI_SOLVER_HPP
