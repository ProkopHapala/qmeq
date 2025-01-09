#pragma once

#include <vector>
#include <cstring>
#include <cmath>
#include "gauss_solver.hpp"
#include "print_utils.hpp"

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
    double* energies;  // State energies [nstates]
    LeadParams* leads; // Lead parameters [nleads]
    double* coupling;  // Coupling matrix elements [nleads * nstates * nstates]
    
    // Constructor
    SolverParams() : energies(nullptr), leads(nullptr), coupling(nullptr) {}
    
    // Copy constructor
    SolverParams(const SolverParams& other) {
        nstates = other.nstates;
        nleads = other.nleads;
        
        // Deep copy arrays
        energies = new double[nstates];
        leads = new LeadParams[nleads];
        coupling = new double[nleads * nstates * nstates];
        
        std::memcpy(energies, other.energies, nstates * sizeof(double));
        std::memcpy(leads, other.leads, nleads * sizeof(LeadParams));
        std::memcpy(coupling, other.coupling, nleads * nstates * nstates * sizeof(double));
    }
    
    // Destructor
    ~SolverParams() {
        delete[] energies;
        delete[] leads;
        delete[] coupling;
    }
    
    // Assignment operator
    SolverParams& operator=(const SolverParams& other) {
        if (this != &other) {
            delete[] energies;
            delete[] leads;
            delete[] coupling;
            
            nstates = other.nstates;
            nleads = other.nleads;
            
            energies = new double[nstates];
            leads = new LeadParams[nleads];
            coupling = new double[nleads * nstates * nstates];
            
            std::memcpy(energies, other.energies, nstates * sizeof(double));
            std::memcpy(leads, other.leads, nleads * sizeof(LeadParams));
            std::memcpy(coupling, other.coupling, nleads * nstates * nstates * sizeof(double));
        }
        return *this;
    }
    
    // Disable move constructor and assignment
    SolverParams(SolverParams&&) = delete;
    SolverParams& operator=(SolverParams&&) = delete;
};

class PauliSolver {
public:
    SolverParams params;
    double* kernel;         // Kernel matrix [nstates * nstates]
    double* rhs;           // Right-hand side vector [nstates]
    double* probabilities; // State probabilities [nstates]
    double* pauli_factors; // Pauli factors for transitions [nleads * nstates * nstates * 2]
    int verbosity;        // Verbosity level for debugging
    std::vector<std::vector<int>> states_by_charge;  // States organized by charge number, like Python's statesdm

    // Count number of electrons in a state
    int count_electrons(int state) {
        return __builtin_popcount(state);
    }

    // Get the site that changed in a transition between two states
    // Returns -1 if more than one site changed or if no site changed
    int get_changed_site(int state1, int state2) {
        int diff = state1 ^ state2;
        if (__builtin_popcount(diff) != 1) {
            return -1;  // More than one site changed or no site changed
        }
        // Find the position of the 1 bit in diff
        return __builtin_ctz(diff);  // Count trailing zeros
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

        if(verbosity > 0) {
            print_vector(states_by_charge, "DEBUG: C++ states_by_charge");
            printf("\nDEBUG: C++ coupling matrix elements in  PauliSolver::init_states_by_charge():\n");
            print_3d_array(params.coupling, params.nleads, n, n, "Lead ");
        }
    }

    // Generate Pauli factors for transitions between states
    void generate_fct() {
        if(verbosity > 0) {
            printf("\nDEBUG: C++ inputs:\n");
            print_vector(params.energies, params.nstates, "State energies (E)");
            
            printf("\nTunneling amplitudes (Tba) (in file pauli_solver.hpp):\n");
            for(int l = 0; l < params.nleads; l++) {
                printf("Lead %d:\n", l);
                print_matrix(&params.coupling[l * params.nstates * params.nstates], 
                           params.nstates, params.nstates);
            }
            
            std::vector<double> mu_vec(params.nleads);
            std::vector<double> temp_vec(params.nleads);
            for(int l = 0; l < params.nleads; l++) {
                mu_vec[l] = params.leads[l].mu;
                temp_vec[l] = params.leads[l].temp;
            }
            print_vector(mu_vec.data(), params.nleads, "Chemical potentials (mu)");
            print_vector(temp_vec.data(), params.nleads, "Temperatures (temp)");
            
            print_vector(states_by_charge, "States by charge");
        }
        
        //exit(0); // DEBUG - We will keep this here until we are sure the leads tunelling amplitudes are correctly pased over the interface

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
                        
                        // Get the site that changed in this transition
                        int changed_site = get_changed_site(c, b);
                        
                        // Calculate coupling strength
                        double coupling_val = params.coupling[l * params.nstates * params.nstates + b * params.nstates + c] * 
                                            params.coupling[l * params.nstates * params.nstates + c * params.nstates + b];
                        
                        // Apply position-dependent coupling for tip (lead 1)
                        if (l == 1 && changed_site > 0) {  // For tip and sites 1,2
                            coupling_val *= 0.09;  // (coeffT * coeffT) = 0.3 * 0.3
                        }
                        
                        // Include lead parameters
                        const LeadParams& lead = params.leads[l];
                        double fermi = fermi_func(energy_diff, lead.mu, lead.temp);
                        
                        // Store factors for both directions
                        // Note: Python's func_pauli multiplies by 2π and coupling already includes gamma/π
                        pauli_factors[idx + 0] = coupling_val * fermi * 2 * PI;         // Forward
                        pauli_factors[idx + 1] = coupling_val * (1.0 - fermi) * 2 * PI; // Backward
                        
                        if(verbosity > 0) printf("DEBUG: generate_fct() l: %d i: %d j: %d E_diff: %.6f coupling: %.6f fermi: %.6f factors:[ %.6f , %.6f ]\n",   l, c, b, energy_diff, coupling_val, fermi, pauli_factors[idx + 0], pauli_factors[idx + 1]);
                    }
                }
            }
        }
        
        // if(verbosity > 0) {
        //     printf("\nDEBUG: Pauli factors:\n");
        //     for(int l = 0; l < params.nleads; l++) {
        //         printf("Lead %d:\n", l);
        //         for(int dir = 0; dir < 2; dir++) {
        //             printf("Direction %d:\n", dir);
        //             print_matrix(&pauli_factors[l * params.nstates * params.nstates * 2 + dir], params.nstates, params.nstates);
        //         }
        //     }
        // }
        
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
        
        // if(verbosity > 0) {
        //     printf("DEBUG: generate_coupling_terms() state:%d diagonal:%.6f\n", state, diagonal);
        //     print_matrix(kernel, n, n, "Phase 2 - After processing state");
        // }
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
        
        // Replace first row with normalization condition
        for(int j = 0; j < n; j++) {
            kernel[j] = 1.0;  // First row is all 1s
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
        
        // Set up RHS vector [1, 0, ..., 0]
        double* rhs = new double[n];
        rhs[0] = 1.0;
        std::fill(rhs + 1, rhs + n, 0.0);
        
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
        printf("DEBUG: PauliSolve() DONE verbosity=%i \n", verbosity);
    }

    ~PauliSolver() {
        delete[] kernel;
        delete[] rhs;
        delete[] probabilities;
        delete[] pauli_factors;
    }

    // Solve the master equation
    void solve() {
        generate_kern();  // First generate the kernel matrix
        solve_kern();     // Then solve it
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
