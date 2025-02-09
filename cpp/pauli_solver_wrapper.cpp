#include "pauli_solver.hpp"
#include <cstdio>
#include "print_utils.hpp"

extern "C" {

// Create a solver instance
void* create_pauli_solver(int nstates, int nleads, 
                         double* energies, double* tunneling_amplitudes,
                         double* lead_mu, double* lead_temp, double* lead_gamma,
                         int verbosity = 0) {
    SolverParams params;
    params.nstates = nstates;
    params.nleads = nleads;
    
    // Allocate and copy energies
    params.energies = new double[nstates];
    std::memcpy(params.energies, energies, nstates * sizeof(double));
    
    // Allocate and set up leads
    params.leads = new LeadParams[nleads];
    for(int i = 0; i < nleads; i++) {
        params.leads[i].mu = lead_mu[i];
        params.leads[i].temp = lead_temp[i];
        params.leads[i].gamma = lead_gamma[i];
    }
    
    // Allocate and copy coupling matrix elements
    params.coupling = new double[nleads * nstates * nstates];
    for(int l = 0; l < nleads; l++) {
        for(int i = 0; i < nstates; i++) {
            for(int j = 0; j < nstates; j++) {
                // Match Python's memory layout
                int idx = (l * nstates + j) * nstates + i;
                params.coupling[l * nstates * nstates + i * nstates + j] = tunneling_amplitudes[idx];
            }
        }
    }
    
    // Debug print tunneling amplitudes
    // if(verbosity > 0) {
    //     printf("DEBUG: tunneling amplitudes after conversion ( in file pauli_solver_wrapper.cpp ):\n");
    //     print_3d_array(tunneling_amplitudes, nleads, nstates, nstates, "Tuneling amplitudes for Lead ");
    // }
    
    PauliSolver* solver = new PauliSolver(params, verbosity);
    
    return solver;
}

// Solve the master equation
void solve_pauli(void* solver_ptr) {
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    solver->solve();
}

// Get the kernel matrix
void get_kernel(void* solver_ptr, double* out_kernel) {
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    const double* kernel = solver->get_kernel();
    int n = solver->params.nstates;
    std::memcpy(out_kernel, kernel, n * n * sizeof(double));
}

// Get the probabilities
void get_probabilities(void* solver_ptr, double* out_probs) {
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    const double* probs = solver->get_probabilities();
    int n = solver->params.nstates;
    std::memcpy(out_probs, probs, n * sizeof(double));
}

// Calculate current through a lead
double calculate_current(void* solver_ptr, int lead_idx) {
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    return solver->generate_current(lead_idx);
}

// Cleanup
void delete_pauli_solver(void* solver_ptr) {
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    delete solver;
}

}
