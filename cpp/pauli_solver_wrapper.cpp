#include "pauli_solver.hpp"
#include <cstdio>

extern "C" {

PauliSolver* create_pauli_solver(
    int nsingle,
    double eps1, double eps2, double eps3,
    double t, double W, double VBias,
    double GammaS, double GammaT,
    double coeffE, double coeffT,
    double muS, double muT, double Temp,
    int verbosity
) {
    // Create system parameters
    SystemParams params(
        nsingle,
        eps1, eps2, eps3,
        t, W, VBias,
        GammaS, GammaT,
        coeffE, coeffT,
        muS, muT, Temp
    );
    
    return new PauliSolver(params, verbosity);
}

void solve_pauli(PauliSolver* solver) {
    solver->solve();
}

const double* get_kernel(PauliSolver* solver) {
    return solver->get_kernel();
}

const double* get_probabilities(PauliSolver* solver) {
    return solver->get_probabilities();
}

const double* get_pauli_factors(PauliSolver* solver) {
    return solver->get_pauli_factors();
}

const double* get_rhs(PauliSolver* solver) {
    return solver->get_rhs();
}

void delete_pauli_solver(PauliSolver* solver) {
    delete solver;
}

}
