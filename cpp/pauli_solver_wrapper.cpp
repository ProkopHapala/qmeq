#include "pauli_solver.hpp"
#include <cstdio>

extern "C" {

// Create a solver instance
void* create_pauli_solver(int nstates, int nleads, 
                         double* energies, double* tunneling_amplitudes,
                         double* lead_mu, double* lead_temp, double* lead_gamma) {
    SystemParams params;
    params.nstates = nstates;
    params.nleads = nleads;
    
    // Allocate and copy energies
    params.energies = new double[nstates];
    std::memcpy(params.energies, energies, nstates * sizeof(double));
    
    // Allocate and copy tunneling amplitudes
    params.tunneling_amplitudes = new double[nleads * nstates * nstates];
    std::memcpy(params.tunneling_amplitudes, tunneling_amplitudes, 
                nleads * nstates * nstates * sizeof(double));
    
    // Allocate and set up leads
    params.leads = new LeadParams[nleads];
    for(int i = 0; i < nleads; i++) {
        params.leads[i].mu = lead_mu[i];
        params.leads[i].temp = lead_temp[i];
        params.leads[i].gamma = lead_gamma[i];
    }
    
    return new PauliSolver(params);
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
    int n = solver->get_nstates();
    std::memcpy(out_kernel, kernel, n * n * sizeof(double));
}

// Get the probabilities
void get_probabilities(void* solver_ptr, double* out_probs) {
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    const double* probs = solver->get_probabilities();
    int n = solver->get_nstates();
    std::memcpy(out_probs, probs, n * sizeof(double));
}

// Calculate current through a lead
double calculate_current(void* solver_ptr, int lead_idx) {
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    return solver->calculate_current(lead_idx);
}

// Cleanup
void delete_pauli_solver(void* solver_ptr) {
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    delete solver;
}

}
