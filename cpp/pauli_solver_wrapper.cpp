#include "pauli_solver.hpp"
#include <cstdio>
#include "print_utils.hpp"

#include <cstdio>  // Make sure this is included

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
        params.leads[i].mu    = lead_mu[i];
        params.leads[i].temp  = lead_temp[i];
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

void* create_pauli_solver_new(int nSingle, int nstates, int nleads, double* Hsingle, double W, double* TLeads, double* lead_mu, double* lead_temp, double* lead_gamma, int* state_order, int verbosity = 0) {
    setvbuf(stdout, NULL, _IONBF, 0);  // Disable buffering for stdout
    //printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    //printf("create_pauli_solver_new() nSingle=%i nstates=%i nleads=%i Hsingle=%p W=%f TLeads=%p lead_mu=%p lead_temp=%p lead_gamma=%p state_order=%p\n", nSingle, nstates, nleads, Hsingle, W, TLeads, lead_mu, lead_temp, lead_gamma, state_order);
    //DEBUG
    SolverParams params;
    params.nSingle = nSingle;
    params.nstates = nstates;
    params.nleads  = nleads;
    params.reallocate(nstates, nleads);
    //DEBUG
    for(int i = 0; i < nleads; i++) {
        params.leads[i].mu    = lead_mu[i];
        params.leads[i].temp  = lead_temp[i];
        params.leads[i].gamma = lead_gamma[i];
    }
    //DEBUG
    for(int i = 0; i < nstates; i++) {
        params.state_order[i] = state_order[i];
    }
    //DEBUG
    //printf("DEBUG: Hsingle = %p\n", Hsingle);
    // Use nSingle for the Hsingle matrix dimensions, not nstates
    params.calculate_state_energies(Hsingle, W);
    //DEBUG
    //printf("DEBUG: params.energies after calculate_state_energies = %p\n", params.energies);
    // Use nSingle for single-particle states
    params.calculate_tunneling_amplitudes(nleads, nstates, nSingle, TLeads);
    //DEBUG
    PauliSolver* solver = new PauliSolver(params, verbosity);
    //DEBUG
    //printf("DEBUG: solver->params.energies after constructor = %p\n", solver->params.energies);
    //printf("DEBUG: solver = %p\n", solver);
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

// Get the energies
void get_energies(void* solver_ptr, double* out_energies) {
    // Print full 64-bit pointer value
    //printf("DEBUG: solver_ptr = %p (full 64-bit address)\n", solver_ptr); 
    // Check if pointer is valid (basic sanity check)    
    PauliSolver* solver = static_cast<PauliSolver*>(solver_ptr);
    //printf("DEBUG: solver = %p\n", solver);
    // Verify solver pointer is valid before using it        
    int n = solver->params.nstates;
    const double* energies = solver->get_energies();
    std::memcpy(out_energies, energies, n * sizeof(double));
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
