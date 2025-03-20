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

### 7. State Ordering

- QmeQ:
  - States ordered by charge first, then by energy
  - Uses `si.statesdm` for charge-based grouping
  - Example order: `[0, 1, 2, 4, 3, 5, 6, 7]`

- C++:
  - Currently uses numerical ordering by binary state value
  - Should implement charge-based sorting similar to QmeQ
  - Requires modifying `SolverParams` constructor to:
    1. Group states by charge
    2. Sort within charge groups by energy
    3. Maintain original QmeQ ordering for matrix compatibility



# What is the problem now ?

OK, I think this problem is somehow related to what I already obsreved before when I was comparing the energies. To match the energies (which are now matching) it was necessary to introduce this state ordering
state_order = [0, 4, 2, 6, 1, 5, 3, 7]
which is basically doing two operation:
1) assuming the states are ordered by charge like this
q=0: [0:'000'],  q=1: [1:'001', 2:'010', 3:'100'],  q=2: [4:'011', 5:'101', 6:'110'], q=3: [7:'111']
shorthand [0 | 1, 2, 3 | 4, 5, 6 | 7]
2) when I considered the bit-reversal I get this [0 | 3, 2, 1 | 6, 5, 4 | 7]
3) But in in ouput from QmeQ this ordering where 3 and 4 is swaped for some resason
Final charge grouping:
Charge 0: states [0] (binary: ['000'])
Charge 1: states [1, 2, 4] (binary: ['001', '010', '100'])
Charge 2: states [3, 5, 6] (binary: ['011', '101', '110'])
Charge 3: states [7] (binary: ['111'])
... not sure why there is this swap in QmeQ
4) so finally I set state ordering to  state_order = [0 | 4, 2, 6 |  1, 5, 3 | 7]
I load it into params.state_order in @pauli_solver.hpp
this makes the energies consistent. But I think it is used only for SolverParams::calculate_state_energies() and it is not used any more for eval_lead_coupling generate_fct (?)
Please check that

then propose how to repait it