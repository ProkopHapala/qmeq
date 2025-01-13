#pragma once

#include <vector>
#include <algorithm>
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
    //int nstates;         // Number of states
    double* kernel;        // Kernel matrix [nstates * nstates]
    double* rhs;           // Right-hand side vector [nstates]
    double* probabilities; // State probabilities [nstates]
    double* pauli_factors; // Pauli factors for transitions [nleads * nstates * nstates * 2]
    double* pauli_factors_compact;  // [nleads][ndm1][2]
    int ndm1;                       // Number of valid transitions (states differing by 1 charge)
    int verbosity;        // Verbosity level for debugging
    std::vector<std::vector<int>> states_by_charge;  // States organized by charge number, like Python's statesdm
    std::vector<int> state_order;      // Maps original index -> ordered index
    std::vector<int> state_order_inv;  // Maps ordered index -> original index

    // from indexing.py of QmeQ
    std::vector<int> lenlst;     // Number of states in each charge sector
    std::vector<int> dictdm;      // State enumeration within charge sectors
    std::vector<int> shiftlst0;  // Cumulative offsets
    std::vector<int> shiftlst1;  // Cumulative offsets
    std::vector<int> mapdm0;      // State pair mapping

    // Count number of electrons in a state
    int count_electrons(int state) {
        return __builtin_popcount(state);
    }

    void count_valid_transitions() {
        ndm1 = 0;
        for(int charge = 0; charge < states_by_charge.size()-1; charge++) {
            ndm1 += states_by_charge[charge+1].size() * states_by_charge[charge].size();
        }
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
    // void init_states_by_charge() {
    //     printf("DEBUG: C++ PauliSolver::init_states_by_charge()\n");
    //     const int n = params.nstates;
    //     int max_charge = 0;
    //     // Find maximum charge number
    //     for(int i = 0; i < n; i++) {
    //         max_charge = std::max(max_charge, count_electrons(i));
    //     }
    //     // Initialize the vector with empty vectors
    //     states_by_charge.resize(max_charge + 1);
    //     // Fill in states for each charge
    //     for(int i = 0; i < n; i++) {
    //         int charge = count_electrons(i);
    //         states_by_charge[charge].push_back(i);
    //     }

    //     if(verbosity > 0) {
    //         printf("\nDEBUG: init_states_by_charge() states_by_charge: \n");
    //         print_vector_of_vectors(states_by_charge );
    //         printf("\nDEBUG: C++ coupling matrix elements in  PauliSolver::init_states_by_charge():\n");
    //         print_3d_array(params.coupling, params.nleads, n, n, "Lead ");
    //     }
    // }

    // Initialize these in init_states_by_charge()
    void init_indexing_maps() {
        const int n = params.nstates;
        lenlst.resize(states_by_charge.size());
        dictdm.resize(n);
        shiftlst0.resize(states_by_charge.size() + 1, 0);
        
        // Fill the mapping arrays following QmeQ's logic
        for(int charge = 0; charge < states_by_charge.size(); charge++) {
            lenlst   [charge]   = states_by_charge[charge].size();
            shiftlst0[charge+1] = shiftlst0[charge] + lenlst[charge] * lenlst[charge];
            
            int counter = 0;
            for(int state : states_by_charge[charge]) {
                dictdm[state] = counter++;
            }
        }

        // Initialize map_dm0
        int total_pairs = shiftlst0.back();  // Total number of possible state pairs
        mapdm0.resize(total_pairs);
        // Fill map_dm0 following QmeQ's logic
        for(int charge = 0; charge < states_by_charge.size(); charge++) {
            for(int b : states_by_charge[charge]) {
                for(int bp : states_by_charge[charge]) {
                    int idx = lenlst[charge] * dictdm[b] + dictdm[bp] + shiftlst0[charge];
                    mapdm0[idx] = idx;  // For now, identity mapping. Add symmetry handling later
                }
            }
        }

        if(verbosity > 0) {
            printf("DEBUG: init_indexing_maps() len_list:    \n"); print_vector(lenlst);
            printf("DEBUG: init_indexing_maps() dict_dm:     \n"); print_vector(dictdm);
            printf("DEBUG: init_indexing_maps() shift_list0: \n"); print_vector(shiftlst0);
            printf("DEBUG: init_indexing_maps() shift_list1: \n"); print_vector(shiftlst1);
            printf("DEBUG: init_indexing_maps() map_dm0:     \n"); print_vector(mapdm0);
        }
    }

    void init_state_ordering() {
        const int n = params.nstates;
        state_order.resize(n);
        state_order_inv.resize(n);
        
        int idx = 0;
        for(int charge = 0; charge < states_by_charge.size(); charge++) {
            for(int state : states_by_charge[charge]) {
                state_order[state] = idx;
                state_order_inv[idx] = state;
                idx++;
            }
        }
        
        if(verbosity > 0) {
            printf("State ordering map: original -> ordered\n");
            for(int i = 0; i < n; i++) {   printf("%d -> %d\n", i, state_order[i]);    }
        }
    }

    void init_states_by_charge() {
        printf("DEBUG: C++ PauliSolver::init_states_by_charge()\n");
        const int n = params.nstates;
        int max_charge = 0;
        
        // First find maximum charge
        for(int i = 0; i < n; i++) {
            max_charge = std::max(max_charge, count_electrons(i));
        }
        states_by_charge.resize(max_charge + 1);
        
        // Fill states in order
        for(int i = 0; i < n; i++) {
            int charge = count_electrons(i);
            states_by_charge[charge].push_back(i);
        }
        
        // Sort states within each charge sector
        for(auto& states : states_by_charge) {
            std::sort(states.begin(), states.end());
        }

        init_state_ordering();
        init_indexing_maps();

        if(verbosity > 0) {
            printf("\nDEBUG: init_states_by_charge() states_by_charge: \n");
            print_vector_of_vectors(states_by_charge);
        }
    }

    void generate_fct_compact() {
        const int n = params.nstates;
        count_valid_transitions();
        pauli_factors_compact = new double[params.nleads * ndm1 * 2]();   // Allocate compact array
        
        // Iterate through charge states
        for(int charge = 0; charge < states_by_charge.size()-1; charge++) {
            int next_charge = charge + 1;
            for(int c : states_by_charge[next_charge]) {
                for(int b : states_by_charge[charge]) {
                    int cb = get_ind_dm1(c, b, charge);  // Get compact index
                    double energy_diff = params.energies[c] - params.energies[b];
                    
                    for(int l = 0; l < params.nleads; l++) {
                        // Calculate coupling
                        double tij = params.coupling[l * n * n + b * n + c];
                        double tji = params.coupling[l * n * n + c * n + b];
                        double coupling_val = tij * tji;
                        
                        // Calculate Fermi factor
                        const LeadParams& lead = params.leads[l];
                        double fermi = fermi_func(energy_diff, lead.mu, lead.temp);
                        
                        // Store in compact format
                        pauli_factors_compact[l * ndm1 * 2 + cb * 2 + 0] = coupling_val * fermi * 2 * PI;
                        pauli_factors_compact[l * ndm1 * 2 + cb * 2 + 1] = coupling_val * (1.0 - fermi) * 2 * PI;
                    }
                }
            }
        }
    }

    // Generate Pauli factors for transitions between states
    void generate_fct() {
        if(verbosity > 0) printf( "\nDEBUG: PauliSolver::%s %s \n", __func__, __FILE__);
        
        const int n  = params.nstates;
        const int n2 = n * n;
        memset(pauli_factors, 0, params.nleads * n * n * 2 * sizeof(double));
        
        // Make sure states are organized by charge
        if(states_by_charge.empty()) { init_states_by_charge();}
        
        if(verbosity > 0) printf( "\nDEBUG: PauliSolver::%s %s \n", __func__, __FILE__);
        // Iterate through charge states (like Python's implementation)
        for(int charge = 0; charge < states_by_charge.size() - 1; charge++) {
            int next_charge = charge + 1;
            
            // Iterate through states in current and next charge state
            for(int c : states_by_charge[next_charge]) {
                for(int b : states_by_charge[charge]) {
                    double energy_diff = params.energies[c] - params.energies[b];
                    
                    // For each lead
                    for(int l = 0; l < params.nleads; l++) {
                        const int idx = l * n2 * 2 + c * n * 2 + b * 2;
                        
                        // Get the site that changed in this transition
                        int changed_site = get_changed_site(c, b);
                        
                        // Calculate coupling strength
                        double tij =  params.coupling[l * n2 + b * n + c];
                        double tji =  params.coupling[l * n2 + c * n + b];
                        double coupling_val = tij * tji; 
                                            
                        // Include lead parameters
                        const LeadParams& lead = params.leads[l];
                        double fermi = fermi_func(energy_diff, lead.mu, lead.temp);
                        
                        // Store factors for both directions
                        // Note: Python's func_pauli multiplies by 2π and coupling already includes gamma/π
                        pauli_factors[idx + 0] = coupling_val * fermi * 2 * PI;         // Forward
                        pauli_factors[idx + 1] = coupling_val * (1.0 - fermi) * 2 * PI; // Backward
                        
                        if(verbosity > 0){
                            //if( (l==1) && (c==3) && (b==1) ){
                            //    printf("DEBUG: generate_fct() l: %d i: %d j: %d E_diff: %.6f coupling: %.6f tij: %.6f tji: %.6f fermi: %.6f factors:[ %.6f , %.6f ]\n",   l, c, b, energy_diff, coupling_val, tij, tji, fermi, pauli_factors[idx + 0], pauli_factors[idx + 1]);
                            //}
                            printf("DEBUG: generate_fct() l: %d i: %d j: %d E_diff: %.6f coupling: %.6f fermi: %.6f factors:[ %.6f , %.6f ]\n",   l, c, b, energy_diff, coupling_val, fermi, pauli_factors[idx + 0], pauli_factors[idx + 1]);
                        }
                    }
                }
            }
        }
        
    }

    int get_ind_dm0(int b, int bp, int charge) {
        //if( maptype==1 ){ return map_dm0[idx]} 
        int idx = lenlst[charge]*dictdm[b] + dictdm[bp] + shiftlst0[charge];
        return mapdm0[ idx ];
        //return len_list[charge] * dict_dm[b] + dict_dm[bp] + shift_list0[charge]; 
    }

    int get_ind_dm1(int c, int b, int bcharge) { // For transitions between states with charge difference of 1
        //int ccharge = bcharge + 1;
        //return len_list[ccharge] * dict_dm[c] + dict_dm[b] + shift_list0[bcharge];
        return lenlst[bcharge]*dictdm[c] + dictdm[b] + shiftlst1[bcharge];
    }

    /// @brief Adds a real value (fctp) to the matrix element connecting the states bb and aa in the Pauli kernel. 
    /// In addition, adds another real value (fctm) to the diagonal element kern[bb, bb].
    /// @param fctm Value to be added to the diagonal element kern[bb, bb].
    /// @param fctp Value to be added to the matrix element connecting states bb and aa.
    /// @param bb Index of the first state.
    /// @param aa Index of the second state.
    /// @note Modifies the internal kernel matrix.
    void set_matrix_element_pauli(double fctm, double fctp, int bb, int aa) {
        int n = params.nstates;
        kernel[bb*n+bb] += fctm; // diagonal
        kernel[bb*n+aa] += fctp; // off-diagonal
    }

    inline int index_paulifct        (int l, int i, int j){ return 2*( j + params.nstates*( i + l*params.nstates )); }
    inline int index_paulifct_compact(int l, int i       ){ return 2*( i + ndm1*l); }

    void generate_coupling_terms(int b) {
        //if(verbosity > 0){ printf("#\n ======== generate_coupling_terms() b: %i \n", b ); }
        const int n = params.nstates;
        int Q = count_electrons(b);
        //const int bb = b * n + b;
        const int bb = b;

        if(verbosity > 0){  printf("\nC++ pauli_solver.hpp ======== generate_coupling_terms() b: %i Q: %i \n", b, Q );  }

        int n2 = n * n;

        if( Q>0 ){ // Handle transitions from lower charge states (a -> b)
            int Qlower=Q-1;
            if(verbosity > 0){ printf("generate_coupling_terms() Q-1 states: " );  print_vector( states_by_charge[Qlower].data(), states_by_charge[Qlower].size()); printf("\n"); }          // for (int a : states_by_charge[Q-1]) printf("%i ", a); printf("\n");
            
            for (int a : states_by_charge[Qlower]) {
                //if (get_changed_site(b, a) == -1) continue;

                int aa = a; // Original
                //int aa = get_ind_dm0(a, a, Qlower);
                //int ba = get_ind_dm1(b, a, Qlower);
                
                double fctm = 0.0, fctp = 0.0;
                for (int l = 0; l < params.nleads; l++) {
                    //int idx = l * n2 * 2 + b * n * 2 + a * 2;
                    int idx = index_paulifct( l, b, a);
                    fctm -= pauli_factors[idx + 1];
                    fctp += pauli_factors[idx + 0];
                }
                //int aa = a * n + a;
                
                if(verbosity > 0){ printf("LOWER [%i,%i] fctm: %.6f fctp: %.6f    bb: %i aa: %i \n", b, a, fctm, fctp, bb, aa); }
                set_matrix_element_pauli(fctm, fctp, bb, aa );
            }
        }        
        if( Q<states_by_charge.size()-1 ){ // Handle transitions to higher charge states (b -> c) 
            int Qhigher=Q+1;
            if(verbosity > 0){ printf("generate_coupling_terms() Q+1 states: " );  print_vector( states_by_charge[Qhigher].data(), states_by_charge[Qhigher].size() ); printf("\n"); } 
            for (int c : states_by_charge[Qhigher]) {
                //if (get_changed_site(b, c) == -1) continue;

                int cc = c;
                //int cc = si.get_ind_dm0(c, c, Qhigher );
                //int cb = si.get_ind_dm1(c, b, Q       );
                double fctm = 0.0, fctp = 0.0;
                for (int l = 0; l < params.nleads; l++) {
                    //int idx = l * n2 * 2 + c * n * 2 + b * 2;
                    int idx = index_paulifct( l, c, b );
                    fctm -= pauli_factors[idx + 0];
                    fctp += pauli_factors[idx + 1];
                }
                //int cc = c * n + c;
                
                if(verbosity > 0){ printf("HIGHER [%i,%i] fctm: %.6f fctp: %.6f    bb: %i aa: %i \n", b, c, fctm, fctp, bb, cc); }
                set_matrix_element_pauli( fctm, fctp, bb, cc );
            }
        }
        //if(verbosity > 0) { printf( "generate_coupling_terms() b: %i kernel: \n", b ); print_matrix(kernel, n, n); }

    }

    void generate_coupling_terms_compact(int b) {
        //if(verbosity > 0){ printf("#\n ======== generate_coupling_terms() b: %i \n", b ); }
        const int n = params.nstates;
        int Q = count_electrons(b);
        //const int bb = b * n + b;
        const int bb = b;

        if(verbosity > 0){  printf("\nC++ pauli_solver.hpp ======== generate_coupling_terms() b: %i Q: %i \n", b, Q );  }

        int n2 = n * n;

        if( Q>0 ){ // Handle transitions from lower charge states (a -> b)
            int Qlower=Q-1;
            if(verbosity > 0){ printf("generate_coupling_terms() Q-1 states: " );  print_vector( states_by_charge[Qlower].data(), states_by_charge[Qlower].size()); printf("\n"); }          // for (int a : states_by_charge[Q-1]) printf("%i ", a); printf("\n");
            
            for (int a : states_by_charge[Qlower]) {
                //if (get_changed_site(b, a) == -1) continue;

                //int aa = a; // Original
                int aa = get_ind_dm0(a, a, Qlower);
                int ba = get_ind_dm1(b, a, Qlower);
                
                double fctm = 0.0, fctp = 0.0;
                for (int l = 0; l < params.nleads; l++) {
                    //int idx = l * n2 * 2 + b * n * 2 + a * 2;
                    int idx = index_paulifct_compact( l, ba );
                    fctm -= pauli_factors_compact[idx + 1];
                    fctp += pauli_factors_compact[idx + 0];
                }
                //int aa = a * n + a;
                
                if(verbosity > 0){ printf("LOWER [%i,%i] fctm: %.6f fctp: %.6f    bb: %i aa: %i \n", b, a, fctm, fctp, bb, aa); }
                set_matrix_element_pauli(fctm, fctp, bb, aa );
            }
        }        
        if( Q<states_by_charge.size()-1 ){ // Handle transitions to higher charge states (b -> c) 
            int Qhigher=Q+1;
            if(verbosity > 0){ printf("generate_coupling_terms() Q+1 states: " );  print_vector( states_by_charge[Qhigher].data(), states_by_charge[Qhigher].size() ); printf("\n"); } 
            for (int c : states_by_charge[Qhigher]) {
                //if (get_changed_site(b, c) == -1) continue;

                //int cc = c;
                int cc = get_ind_dm0(c, c, Qhigher );
                int cb = get_ind_dm1(c, b, Q       );
                double fctm = 0.0, fctp = 0.0;
                for (int l = 0; l < params.nleads; l++) {
                    //int idx = l * n2 * 2 + c * n * 2 + b * 2;
                    int idx = index_paulifct_compact( l, cb );
                    fctm -= pauli_factors_compact[idx + 0];
                    fctp += pauli_factors_compact[idx + 1];
                }
                //int cc = c * n + c;
                
                if(verbosity > 0){ printf("HIGHER [%i,%i] fctm: %.6f fctp: %.6f    bb: %i aa: %i \n", b, c, fctm, fctp, bb, cc); }
                set_matrix_element_pauli( fctm, fctp, bb, cc );
            }
        }
        //if(verbosity > 0) { printf( "generate_coupling_terms() b: %i kernel: \n", b ); print_matrix(kernel, n, n); }

    }

    void normalize_kernel() {
        const int n = params.nstates;
        // Set first row to all ones (like Python)
        for(int j = 0; j < n; j++) { kernel[j] = 1.0; }        
        // if(verbosity > 0) {
        //     printf("Phase 2 - After normalization\n");
        //     print_matrix(kernel, n, n, nullptr, "%.6g");
        // }
    }

    // Generate kernel matrix
    void generate_kern() {
        if(verbosity > 0) printf("\nDEBUG: generate_kern() Building kernel matrix...\n");
        const int n = params.nstates;
        generate_fct();
        std::fill(kernel, kernel + n * n, 0.0);

        for(int state = 0; state < n; state++) { 
            int b = state_order_inv[state];
            generate_coupling_terms(b); 
        }
        if(verbosity > 0) { printf( "generate_kern(): kernel: \n" ); print_matrix(kernel, n, n); }
        //normalize_kernel();
        //if(verbosity > 0) { print_matrix(kernel, n, n, "Phase 2 - After normalization"); }
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
