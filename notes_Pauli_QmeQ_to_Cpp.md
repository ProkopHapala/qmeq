# Notes on Porting QmeQ Pauli Solver to C++

## Task Overview
Our goal is to make the C++ solver in `cpp/pauli_solver.hpp` reproduce the results of the Python QmeQ solver, which is mostly implemented in `qmeq/approach/base/pauli.py` but also in related classes like `qmeq/indexing.py` and `qmeq/qdot.py`.

## Current Progress
Energies of 8 = 2^3 many-body states are calculated properly. Now we need to investigate the building of the kernel matrix with coupling between these 8 states. Debug outputs are printed with verbosity > 3.

## Relevant Functions and Files
- `compare_solvers.py`: Script to compare the solvers
- `compare_solvers.log`: Log file with debug outputs
- C++
    - `cpp/pauli_solver.hpp`: C++ implementation of the solver
    - `cpp/gauss_solver.hpp`: C++ implementation of the least squares solver
    - `cpp/pauli_solver_wrapper.cpp`: C++ wrapper for the solver
    - `pauli_solver_lib.py`: Python binding which loads the C++ library compiled from `cpp/pauli_solver_wrapper.cpp` using ctypes (see also `compare_solvers.py`)
- Python QmeQ
    - `qmeq/approach/base/pauli.py`: Main Python implementation of the solver
    - `qmeq/indexing.py`: Indexing-related functions
    - `qmeq/qdot.py`: Quantum dot-related functions
    - `qmeq/leadstun.py`: creates Tunneling amplitudes (Tba) interaction of many-body states with leads

## Dimension Variables
- `nSingle`: Number of single-particle states
- `nMany`: Number of many-body states (2^nSingle)
- `nLeads`: Number of leads
- `nCharge`: Number of charge sectors (nSingle + 1)

## Pauli Calculation Scheme in QmeQ

### 1. State Indexing and Ordering
- QmeQ:
  - `StateIndexing` class (`qmeq/indexing.py`) manages state ordering
  - `construct_chargelst()` (`qmeq/indexing.py`) groups states by charge
  - the ordering seems fishy. It seems that states are ordered by charge like this: `[ [0], [1, 2, 4], [3, 5, 6], [7] ]` (i.e. `[['000'], ['001', '010', '100'], ['011', '101', '110'], ['111'] ]`) 
     - but if we use this ordering in C++ we get results for energies that are different from QmeQ (reshuffled states)
- C++ equivalent: implicit in state ordering
    - in C++ we can read external state ordering from `state_order` when calling `create_pauli_solver_new()` in (see `compare_solvers.py`)
    - currently is seems that state_order = [0, 4, 2, 6, 1, 5, 3, 7] reproduces energies of QmeQ

### 2. Energy Calculation

- QmeQ:
  - `construct_Ea_manybody()` (`qmeq/qdot.py`) calculates many-body energies
  - as input use Single-particle Hamiltonian `Hsingle` shape: [nSingle, nSingle]
  - and Coulomb interaction `W` between occupied sites
- C++:
  - `calculate_state_energy()` (`pauli_solver.hpp`) is equivalent of QmeQ `construct_Ea_manybody()`

### 3. Lead Coupling

- QmeQ:
- input is Single-particle tunneling amplitudes `tleads` shape:[nLeads ,nSingle]
- Many-body hopping amplitudes between many-body states due to interaction with leads are stored in`Tba` shape:[nLeads,nMany,nMany]
  - Created by `construct_Tba()` in `LeadsTunneling` class (`qmeq/leadstun.py`)
  - Includes fermion sign factors for electron addition/removal (+/-1)
- C++:
  - `SolverParams::calculate_tunneling_amplitudes()` and `SolverParams::eval_lead_coupling()` (in `pauli_solver.hpp`) is equivalent of QmeQ `construct_Tba()`

### 4. Kernel Matrix Generation

- QmeQ:
  - Rates between states stored in matrix calculated in `ApproachPauli.generate_fct()` (`qmeq/approach/base/pauli.py`)
- Kernel matrix `K` shape:[nMany,nMany] built in `ApproachPauli.generate_kern()`
  - Depends on:
    - State energies `E` [nMany]
    - Tunneling amplitudes `Tba` shape:[nLeads,nMany,nMany]
    - Fermi functions of leads
- C++:
  - `generate_kern()` (`pauli_solver.hpp`) is equivalent of QmeQ `ApproachPauli.generate_kern()`

### 5. State Probabilities Calculation

- QmeQ:
  - Probabilities calculated from steady-state solution of Pauli master equation
  - Kernel matrix `K` modified for probability conservation in `ApproachPauli.generate_kern()` (`qmeq/approach/base/pauli.py`)
    - One row replaced with ones to enforce probability conservation
    - Matrix equation `K·p = 0` solved for steady-state probabilities `p`
- C++:
  - Probabilities calculated similarly in `solve_steady_state()` (`pauli_solver.hpp`)
  - Kernel matrix modification handled in `generate_kern()`
  - Uses same probability conservation approach as QmeQ

### 6. Current Calculation

- QmeQ:
  - Currents calculated using kernel matrix elements
  - Depends on state probabilities and transition rates
  - Calculated in `ApproachPauli.generate_current()` (`qmeq/approach/base/pauli.py`)
- C++:
  - `generate_current()` (`pauli_solver.hpp`) is equivalent of QmeQ `ApproachPauli.generate_current()`

