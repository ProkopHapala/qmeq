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
