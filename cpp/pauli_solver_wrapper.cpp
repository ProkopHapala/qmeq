#include "pauli_solver.hpp"
#include "print_utils.hpp"

extern "C" {

PauliSolver* create_pauli_solver(
    int nstates, int nsingle,
    const double* energies,
    double eps1, double eps2, double eps3,
    double t, double W, double VBias,
    double GammaS, double GammaT,
    double coeffE, double coeffT,
    const double* lead_mu,
    const double* lead_temp,
    int verbosity
) {
    // Create parameters struct
    SystemParams params;
    params.nstates = nstates;
    params.nsingle = nsingle;
    params.nleads = 2;  // Fixed for this problem
    
    // Copy energies
    params.energies.resize(nstates);
    std::copy(energies, energies + nstates, params.energies.begin());
    
    // Set single-particle parameters
    params.eps1 = eps1;
    params.eps2 = eps2;
    params.eps3 = eps3;
    params.t = t;
    params.W = W;
    params.VBias = VBias;
    params.GammaS = GammaS;
    params.GammaT = GammaT;
    params.coeffE = coeffE;
    params.coeffT = coeffT;
    
    // Calculate tunneling amplitudes VS and VT
    params.VS = sqrt(params.GammaS * params.W / (2 * 3.14159265358979323846));
    params.VT = sqrt(params.GammaT * params.W / (2 * 3.14159265358979323846));
    
    // Copy lead parameters
    params.leads.resize(params.nleads);
    for(int l = 0; l < params.nleads; l++) {
        params.leads[l].mu = lead_mu[l];
        params.leads[l].temp = lead_temp[l];
    }
    
    return new PauliSolver(params, verbosity);
}

void delete_pauli_solver(PauliSolver* solver) {
    delete solver;
}

void solve_master_equation(PauliSolver* solver) {
    solver->solve();
}

const double* get_probabilities(PauliSolver* solver) {
    return solver->get_probabilities();
}

const double* get_kernel(PauliSolver* solver) {
    return solver->get_kernel();
}

const double* get_rhs(PauliSolver* solver) {
    return solver->get_rhs();
}

const double* get_pauli_factors(PauliSolver* solver) {
    return solver->get_pauli_factors();
}

}
