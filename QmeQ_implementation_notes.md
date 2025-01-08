# QmeQ Implementation Notes

## State Handling and Tunneling Amplitudes

### State Indexing System
- QmeQ uses a sophisticated state indexing system (`StateIndexing` class) that organizes states by charge sectors
- States are grouped in `statesdm` by charge number
- This organization is crucial for tracking allowed transitions between states
- The indexing system supports both Fock basis and eigenstate basis

### Tunneling Amplitudes Calculation
Key components in `leadstun.py`:

1. **Single-particle Tunneling**
   - Uses `tleads` dictionary: `tleads[(lead, state)] = tunneling_amplitude`
   - Allows flexible coupling configuration for each lead-state pair

2. **Many-body Tunneling Matrix**
   - Built using `construct_Tba` function
   - Matrix elements: `Tba[lead, final_state, initial_state]`
   - Handles both electron addition and removal:
     ```python
     # Adding electron
     Tba[lead, ind, j1] += fsign * tamp
     # Removing electron
     Tba[lead, ind, j1] += fsign * np.conj(tamp)
     ```
   - Fermion sign calculation: `fsign = np.power(-1, sum(state[0:j2]))`

3. **Basis Rotation**
   - Initial tunneling amplitudes are calculated in Fock basis
   - Can be rotated to eigenstate basis using `rotate_Tba` function
   - Rotation uses eigenvectors from quantum dot Hamiltonian diagonalization

### Important Functions

1. `construct_Tba` (leadstun.py)
   - Builds many-body tunneling amplitude matrix from single-particle amplitudes
   - Handles proper sign conventions and state transitions

2. `rotate_Tba` (leadstun.py)
   - Rotates tunneling amplitudes from Fock basis to eigenstate basis
   - Uses eigenvectors from system Hamiltonian

3. `generate_fct` (pauli.py)
   - Calculates Pauli master equation factors
   - Uses energy differences and tunneling amplitudes
   - Applies Fermi functions for lead couplings

### Critical Implementation Details

1. **Matrix Organization**
   - Consistent indexing: `[lead, final_state, initial_state]`
   - Proper handling of Hermitian conjugation through complex conjugation
   - No explicit negative signs for conjugate terms

2. **Sign Conventions**
   - Fermion signs based on occupied states before the transition site
   - Complex conjugation for electron removal processes
   - Signs must match between tunneling amplitudes and state transitions

3. **State Transitions**
   - Only allowed between states differing by one electron
   - Proper tracking of charge sectors for valid transitions
   - Energy conservation in transition rates

## Coulomb Interaction Handling

### Implementation Details

The Coulomb interaction in QmeQ is implemented as a many-body operator in the form:

```
H_coulomb = U * c†_m c†_n c_k c_l
```

where:
- U is the interaction strength
- c†, c are creation/annihilation operators
- m,n,k,l are single-particle state indices

Key implementation points:

1. **Matrix Construction** (`qmeq/qdot.py:construct_ham_coulomb`):
   - Operates on many-body Fock states, not single-particle states
   - Constructs a many-body Hamiltonian matrix of dimension (N_many x N_many)
   - For charge=1 sector, this is a 3x3 matrix
   - For full space, this would be an 8x8 matrix (2^3 states)

2. **Matrix Elements**:
   - Only connects Fock states that differ by exactly two electrons being moved
   - Preserves total electron number
   - Includes proper fermionic signs for electron reordering

3. **Diagonalization**:
   - The resulting Coulomb Hamiltonian is diagonal in the Fock basis
   - No quantum superpositions are created between different sites
   - The Coulomb interaction only affects state energies based on occupation numbers

### Important Implementation Notes

1. **Hamiltonian Structure**:
   - The many-body Hamiltonian is constructed separately for each charge sector
   - For charge=1 sector with 3 single-particle states, we see a 3x3 matrix
   - The matrix is diagonal because Coulomb interaction preserves occupation numbers
   - Example from debug output:
     ```
     DEBUG: construct_manybody_eigenstates() ham:
     [[-10. +0.j   0. +0.j   0. +0.j]
      [  0. +0.j -10. +0.j   0. +0.j]
      [  0. +0.j   0. +0.j -33.6+0.j]]
     ```

2. **Diagonalization Process**:
   - While `diagonalise()` is called, it's effectively a no-op when there's no hopping
   - The eigenvectors are trivial (identity matrix) because:
     - No hopping terms between sites (preserves locality)
     - Coulomb interaction only shifts energies based on occupations
     - No mechanism to create superpositions between different Fock states

3. **Physical Interpretation**:
   - Coulomb interaction U affects energies but not wavefunctions
   - No bonding/anti-bonding states are created
   - States remain localized because interaction depends only on occupation numbers
   - Phase relationships between sites are not affected

The implementation can be found in:
- Main function: `construct_ham_coulomb()` in `qmeq/qdot.py`
- Called by: `construct_manybody_eigenstates()` in same file
- Debug output shows the diagonal nature of the resulting matrix

## System Setup and Initialization

1. **Builder Class**
   - Initializes the full system setup
   - Handles state indexing, quantum dot, and leads
   - Creates appropriate approach (Pauli, Lindblad, etc.)

2. **Quantum Dot**
   - Stores single-particle and interaction Hamiltonians
   - Handles state energies and eigenvectors
   - Manages basis transformations

3. **Leads**
   - Stores lead parameters (chemical potentials, temperatures)
   - Manages tunneling couplings
   - Handles lead-system interactions

## Important Variables

1. **System Parameters**
   - `nsingle`: Number of single-particle states
   - `nleads`: Number of leads
   - `statesdm`: States organized by charge sectors

2. **Tunneling Parameters**
   - `tleads`: Single-particle tunneling amplitudes
   - `Tba`: Many-body tunneling amplitude matrix
   - `mulst`, `tlst`, `dlst`: Lead chemical potentials, temperatures, bandwidths

3. **Quantum Dot Parameters**
   - `hsingle`: Single-particle Hamiltonian
   - `coulomb`: Two-particle interaction terms
   - `Ea`: State energies
