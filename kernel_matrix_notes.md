# Kernel Matrix Calculation Notes

## Task Overview
Our goal is to make the C++ solver in `cpp/pauli_solver.hpp` reproduce the results of the Python QmeQ solver, which is mostly implemented in `qmeq/approach/base/pauli.py` but also in related classes like `qmeq/indexing.py` and `qmeq/qdot.py`.

## Current Progress
Energies of 8 = 2^3 many-body states are calculated properly. Now we need to investigate the building of the kernel matrix with coupling between these 8 states. Debug outputs are printed with verbosity > 3.

## Relevant Functions and Files
- `cpp/pauli_solver.hpp`: C++ implementation of the solver
- `qmeq/approach/base/pauli.py`: Main Python implementation of the solver
- `qmeq/indexing.py`: Indexing-related functions
- `qmeq/qdot.py`: Quantum dot-related functions
- `compare_solvers.py`: Script to compare the solvers
- `compare_solvers.log`: Log file with debug outputs

## Construct of lead tunneling amplitudes

Builds coupling matrix between leads and transitions between many-body states

- In QmeQ Python it is called `leads.Tba` of size [nleads,nmany,nmany] 
   - it is created in `qmeq/leadstun.py`, mainly `construct_Tba` function

- In C++ it is called `SolverParams::coupling` and it is created in SolverParams::eval_lead_coupling()

## Detailed Discrepancy Analysis

Let's analyze key discrepancies from the compare_solvers.log:

1. **Lead 1 Coupling Factors** (Critical Issue):
```log
# QmeQ (Python)
DEBUG: generate_fct() l:1 i:4 j:0 coupling:0.015915
DEBUG: generate_fct() l:1 i:3 j:1 coupling:0.001432

# C++ 
DEBUG: generate_fct() l:1 i:4 j:0 coupling:0.000000
DEBUG: generate_fct() l:1 i:3 j:1 coupling:0.000000
```
**Root Cause**: Lead 1 coupling parameters not being properly passed to C++ solver

2. **Kernel Matrix Differences**:
```log
# QmeQ Python kernel:
[[-1.318       0.409       0.409       0.5        ...]
 [ 0.409      -0.909       0.          0.         ...]]

# C++ kernel: 
[[-0.899       0.284       0.284       0.355      ...]
 [ 0.284      -0.625       0.          0.         ...]]
```
**Pattern**: C++ values ≈ 0.7× Python values → Missing frequency dependence?

3. **State Transitions Not Registered**:
```log
# QmeQ
DEBUG: generate_coupling_terms() state:1 other:3 rate:0.400000

# C++
DEBUG: generate_coupling_terms() state:1 other:3 rate:0.000000
```
**Possible Cause**: Incorrect charge state mapping in states_by_charge



**Step-by-Step Correction Plan**:

1. **Verify Coupling Parameter Loading**:
```cpp
// In PauliSolver initialization:
printf("DEBUG: Lead 1 coupling matrix:\n");
for(int i=0; i<nstates; i++){
    for(int j=0; j<nstates; j++){
        printf("%.6f ", params.coupling[1*nstates*nstates + i*nstates + j]);
    }
    printf("\n");
}
```

2. **Add Transition Energy Debug**:
```cpp
// In generate_fct():
printf("DEBUG: Transition %d->%d E_diff=%.6f Tba=%.6f Tba_rev=%.6f\n",
    b, c, energy_diff, tij, tji);
```

3. **Validate Charge State Groups**:
```cpp
// In init_states_by_charge():
printf("DEBUG: Charge group %d states: [", charge);
for(int s : states_by_charge[charge]) printf("%d ", s);
printf("]\n");
```

4. **Add Kernel Assembly Checks**:
```cpp
// In generate_kern():
printf("DEBUG: Assembling kernel row %d (state %d)\n", bb, b);
printf("DEBUG: Row values: ");
for(int i=0; i<nstates; i++) printf("%.6f ", kern[bb*nstates + i]);
printf("\n");
