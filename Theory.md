
# Formulating the Pauli Master Equation for a Triangular Molecular System

This document provides a "cookbook" guide to formulating the Pauli Master Equation (PME) for a system of three molecules arranged in a triangle, coupled to a substrate (lead 1) and an STM tip (lead 2). It's designed to help you implement your own PME solver in C/C++ or Julia, based on the concepts used in the QmeQ package.

## System Description

We consider three molecules forming a symmetric triangular ring on a substrate. Each molecule is strongly coupled to the substrate and weakly coupled to an STM tip.

*   **On-site Coulomb Repulsion:** Each molecule can hold only one electron due to strong on-site Coulomb repulsion.
*   **Basis States:** The system's basis states can be represented as |abc⟩, where a, b, and c are either 0 (unoccupied) or 1 (occupied), representing the occupation of each molecule. Considering the strong on-site Coulomb repulsion, each molecule can be at max singly occupied. We have therefore $2^3=8$ relevant states: |000⟩, |100⟩, |010⟩, |001⟩, |011⟩, |101⟩, |110⟩, |111⟩.
*   **Transport:** Electrons tunnel between the STM tip and the molecules.

## Goal

Our goal is to derive the Pauli Master Equation (PME) for this system, solve it to find the steady-state probabilities of each basis state, and then calculate the current through the system as a function of the applied voltage.

## Pauli Master Equation (PME) Formulation

The PME describes the time evolution of the probabilities P{i} of finding the system in each basis state |i⟩. It's given by:

$$
dp_i/dt = \sum [W_{ji} p_j - W_{ij}p_i]
$$

where:

*   $p_i$ is the probability of the system being in state |i⟩.
*   $W_{ji}$ is the transition rate from state |j⟩ to state |i⟩.

# Steps to Formulate the PME

## Step 1: Define the Hamiltonian

The total Hamiltonian of the system is:

$$
H = H_{dot} + H_{leads} + H_{tunneling}
$$

###  (Quantum Dot Hamiltonian) $H_{dot}$:
    
$H_{dot}     = H_{single} + H_{coulomb}$ 
$H_{single}  = \sum ε_i d_i⁺ d_i + \sum Ω_{ij} (d_i⁺d_j + d_j⁺d_i)$
$H_{coulomb} = \sum U_{ij} n_in_j $
    
where:
*   $ε_i$ is the on-site energy of molecule i.
*   $d_i⁺$, d{i} are the creation and annihilation operators for an electron on molecule i.
*   $Ω_{ij}$ is the coupling between molecules i and j.
*   $U_{ij}$ is the Coulomb interaction between electrons on molecules i and j.
*   $n_i = d_i⁺d_i$ is the number operator for molecule i.

In QmeQ, `hsingle` (dictionary) stores ε{i} and Ω{ij}, and `coulomb` (dictionary) stores U{ij}.

### Leads Hamiltonian $H_{leads}$ :
    
$H_{leads} = \sum ε_{kα}c_{kα}⁺c_{kα}$
where:
*   $ε_{kα}$ is the energy of an electron in lead α with momentum k.
*   $c_{kα}⁺, c_{kα}$ are the creation and annihilation operators for an electron in lead α with momentum k.

###   Tunneling Hamiltonian $H_{tunneling}$ :
    
$H_{tunneling} = \sum { t_{αi} d_{i}⁺ c_{kα} + H.c. } $

where:
*   $t_{αi}$ is the tunneling amplitude between molecule i and lead α.
*   H.c. denotes the Hermitian conjugate.

In QmeQ, `tleads` (dictionary) stores $t_{αi}$.

## Step 2: Calculate Transition Rates $W_{ji}$

The transition rates $W_{ji}$ are calculated using Fermi's Golden Rule:

$ W_{ji} = (2π/ħ) |⟨i|H_{tunneling}|j⟩|^2 ρ(E_i) f(E_j - E_i) $

where:

*   $|i⟩$ and $|j⟩$ are the final and initial many-body states, respectively.
*   $ρ(E_i)$ is the density of states in the lead at the final state energy.
*   $f(E)$ is the Fermi-Dirac distribution function: $f(E) = 1 / (\exp((E - μ)/kT) + 1)$.
    *   μ is the chemical potential of the lead.
    *   k is the Boltzmann constant.
    *   T is the temperature.

### Simplifications for PME:

1. **Wide-band limit:** Assume a constant density of states ρ(E) in the leads.
2. **Neglect coherences:** In the PME, we neglect coherences between different many-body states, meaning we only consider the diagonal elements of the density matrix (the probabilities).

### Transition Rate Calculation in QmeQ:

QmeQ calculates the transition rates by first expressing H{tunneling} in the many-body eigenbasis of H{dot}:

$H_{tunneling} = \sum T_{ba,kα}|b⟩⟨a|c_{kα} + H.c.$

where |a⟩ and |b⟩ are many-body eigenstates of H{dot}, and T{ba,kα} are the many-body tunneling amplitudes.

Then, the transition rates are calculated using a modified Fermi's Golden Rule, taking into account the different possible transitions between many-body states and the Fermi-Dirac distribution in the leads. The rates are stored in the `system.kern` matrix in QmeQ.

## Step 3: Construct the PME Matrix

For our three-molecule system with one electron per molecule, we have four relevant basis states:

 * $|000⟩$ - all molecules unoccupied - state 0
 * $|100⟩$ - molecule 1 occupied - state 1
 * $|010⟩$ - molecule 2 occupied - state 2
 * $|001⟩$ - molecule 3 occupied - state 3
 * $|110⟩$ - molecules 1 and 2 occupied - state 4
 * $|101⟩$ - molecules 1 and 3 occupied - state 5
 * $|011⟩$ - molecules 2 and 3 occupied - state 6
 * $|111⟩$ - all molecules occupied - state 7

The PME matrix will be an 8x8 matrix. Each element (i, j) of the matrix represents the transition rate from state $|j⟩$ to state $|i⟩$.


### Example: Calculating Matrix Elements

Let's consider a few examples of how to calculate the matrix elements:

*   **$W_{10}$ (Transition from $|000⟩$ to $|100⟩$):** This corresponds to an electron tunneling from the substrate or tip to molecule 1. The rate will depend on the tunneling amplitude t{α1}, the density of states in the lead, and the Fermi-Dirac function.
*   **$W_{01}$ (Transition from $|100⟩$ to $|000⟩$):** This corresponds to an electron tunneling from molecule 1 to the substrate or tip.
*   **$W_{21}$ (Transition from $|100⟩$ to $|010⟩$):** This corresponds to an electron hopping from molecule 1 to molecule 2. This rate will depend on the coupling strength Ω{12} between the molecules.

### Important Considerations for the Triangular System:

*   **Symmetry:** The symmetry of your system will affect the transition rates. For example, if molecules 1 and 2 are symmetrically coupled to the tip, then W{10} (substrate to molecule 1) might be different from W{30} (substrate to molecule 3) if molecule 3 is further away.
*   **"Dark State":** The formation of a "dark state" (or a state with similar properties) is crucial for observing NDR. This state will have a very weak coupling to the tip, effectively blocking transport. You'll need to carefully choose the parameters (on-site energies, inter-molecular couplings, and tunneling amplitudes) to create such a state.

## Step 4: Solve the PME

In the steady state, $dp_i/dt = 0$. So, we need to solve the following system of linear equations:

$ d p_i/dt =  0 = \sum_j [W_{ji}p_j - W_{ij}p_i] $

We also have the normalization condition: 
$ \sum_i p_i = 1 $

These equations can be written in matrix form:
$L p = b$

where:

*   $L$ is the PME matrix (constructed in Step 3). It is represented by `system.kern` in QmeQ.
*   $p$ is a vector containing the probabilities $p_i$. It is represented by `system.phi0` in QmeQ.
*   $b$ is a vector that ensures normalization (usually containing all zeros except for one element equal to 1).

You can use standard linear algebra libraries in C/C++ or Julia to solve this matrix equation and obtain the steady-state probabilities.


### Construct the linear system $L p = b$

1. **Construct the matrix L `system.kern` :** The elements of the matrix L are derived from the transition rates $W_{ji}$. Specifically:
    *  **Off-diagonal elements (i ≠ j):** $L_{ij} = W_{ji}$ (The transition rate from state |j⟩ to state |i⟩).
    *  **Diagonal elements (i = j):** $L_{ii} = - \sum_{j≠i} W_{ij}$ (The negative sum of transition rates *out* of state |i⟩).
2. **Construct the vector p:** This vector simply contains the probabilities $P_i$ of each state |i⟩: 
  $p = [P₀, P₁, P₂, P₃, P₄, P₅, P₆, P₇]$  (for our 8-state system).
3. **Construct the vector b:** This vector is related to the normalization condition ($\sum_i P_i = 1$). In most cases, it's a vector of zeros, except for one element that is set to 1. This element corresponds to an equation that is replaced by the normalization condition.

**Building the Matrix L and Vector b in Code**

Here's a Python code example demonstrating how to build the matrix L and vector b, assuming you have already calculated the transition rates W_{ji} and stored them in a suitable data structure (e.g., a 2D array or a dictionary):

```python
import numpy as np

def build_PME_matrix(W, num_states):
    L = np.zeros((num_states, num_states))
    b = np.zeros((num_states, 1))
    for i in range(num_states):
        for j in range(num_states):
            if i != j:
            L[i, i] -= W[i, j]  # sum of transition rate out from |i> to |j>
    L[-1, :] = 1.0   # Replace the last row of L with the normalization condition (sum of probabilities = 1)
    b[-1   ] = 1.0

    return L, b
```

**Explanation:**

1. **Initialization:** We initialize the `L` matrix and `b` vector with zeros.
2. **Populating `L`:**
    *   The outer loop iterates over rows (i) of the matrix, and the inner loop iterates over columns (j).
    *   **Off-diagonal elements:** `L[i, j]` is assigned the transition rate `W[j, i]` (from state j to state i). We use `.get((j, i), 0)` method for dictionaries in case some rates are not present (then they are assumed to be 0).
    *   **Diagonal elements:** `L[i, i]` is updated by subtracting `W[i, j]` (from state i to state j) in each iteration of the inner loop.
3. **Normalization:** The last row of `L` is set to all ones, and the last element of `b` is set to 1. This replaces the last equation in the system with the normalization condition.
4. **Solving:** Finally, `np.linalg.solve(L, b)` solves the linear system to find the steady-state probabilities `P`.

**Important Considerations:**

*   **Transition Rate Calculation:** The code assumes that you have already calculated the transition rates `W_rates` based on your system's Hamiltonian, tunneling amplitudes, and the Fermi-Dirac distribution. This is the most physics-involved part of the problem.
*   **Data Structure for `W`:** You can use either a 2D NumPy array or a dictionary to store the transition rates. The code example shows how to handle both cases.
*   **Numerical Stability:** Be mindful of numerical stability when solving the linear system, especially for larger systems. You might need to use more sophisticated linear algebra solvers or techniques if you encounter issues.

This detailed explanation and code example should clarify how the PME is formulated as a matrix equation and how you can build the necessary matrices in your C/C++ or Julia implementation. Remember that the accurate calculation of the transition rates is crucial for obtaining meaningful results.




## Step 5: Calculate the Current

The current flowing from lead α can be calculated as:

$ I_α = e \sum [W_{ji}^α p_j - W_{ij}^α P_i ] $

where:

*   $e$ is the electron charge.
*   $W_{ji}^α$ is the transition rate from state $|j⟩$ to state $|i⟩$ involving an electron tunneling from/to lead $α$.

In QmeQ, the current is calculated and stored in `system.current`.






# Code Snippets (Illustrative)

Here are some simplified Python code snippets (similar to what's done in QmeQ) to illustrate the key steps:

```python
import numpy as np

# Define system parameters (example)
e0, e1, e2 = 0.0, 0.1, 0.2  # On-site energies
omega12, omega23, omega31 = 0.05, 0.02, 0.05 # Inter-molecular couplings
U12, U23, U31 = 0.0, 0.0, 0.0 # in our case we neglect it, because we have max 1 electron per molecule
t_substrate1, t_tip1 = 0.1, 0.01  # Tunneling amplitudes to molecule 1
t_substrate2, t_tip2 = 0.1, 0.01  # Tunneling amplitudes to molecule 2
t_substrate3, t_tip3 = 0.1, 0.001  # Tunneling amplitudes to molecule 3, notice that it is smaller than for molecules 1 and 2
mu_substrate, mu_tip = 0.0, 0.5  # Chemical potentials
temp = 0.1  # Temperature

# Construct the Hamiltonian (hsingle and coulomb dictionaries)
hsingle = {(0, 0): e0, (1, 1): e1, (2,2): e2,
           (0, 1): omega12, (1, 0): omega12,
           (1, 2): omega23, (2, 1): omega23,
           (0, 2): omega31, (2, 0): omega31}
coulomb = {} # {(0,1,1,0): U12, (1,2,2,1): U23, (0,2,2,0): U31} we neglect it

# Define lead properties (mulst and tlst dictionaries)
nleads = 2
mulst = {0: mu_substrate, 1: mu_tip}
tlst = {0: temp, 1: temp}

# Define tunneling amplitudes (tleads dictionary)
tleads = {(0, 0): t_substrate1, (1, 0): t_tip1,
          (0, 1): t_substrate2, (1, 1): t_tip2,
          (0, 2): t_substrate3, (1, 2): t_tip3}

# Construct the PME matrix (system.kern)
# ... (This part involves calculating the transition rates Wji
#     based on the Hamiltonian, tunneling amplitudes, and Fermi-Dirac
#     distribution. It's the most involved part of the code and
#     requires careful consideration of the energy differences and
#     selection rules.)

# Solve the PME (L * P = B)
# ... (Use a linear algebra solver to find the steady-state probabilities P)

# Calculate the current (system.current)
# ... (Use the formula for I_alpha and the calculated probabilities)
```

**Important Notes:**

*   The provided code snippets are simplified and for illustrative purposes only. The actual implementation in QmeQ is more complex and handles various optimizations and edge cases.
*   The calculation of the transition rates `Wji` is the most crucial and complex part of the implementation. You'll need to carefully consider the energy differences between the states, the Fermi-Dirac distribution, and the tunneling amplitudes.
*   The specific form of the PME matrix and the current calculation will depend on the details of your system and the approximations you make.
*   Pay close attention to the indexing of states and leads in your implementation.

