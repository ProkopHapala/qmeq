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

### Quantum Dot Hamiltonian $H_{dot}$:
    
$H_{dot}     = H_{single} + H_{coulomb}$ 
$H_{single}  = \sum ε_i d_i⁺ d_i + \sum Ω_{ij} (d_i⁺d_j + d_j⁺d_i)$
$H_{coulomb} = \sum U_{ijkl} d_i⁺d_j⁺d_kd_l$

The Hamiltonian is constructed in two steps (see `qmeq/qdot.py`):
1. `construct_ham_hopping()`: Builds $H_{single}$ in the Fock basis
2. `construct_ham_coulomb()`: Adds the Coulomb interaction terms

The eigenstates |i⟩ are then found by diagonalizing H_{dot}.

### Leads Hamiltonian $H_{leads}$ :
    
$H_{leads} = \sum ε_{kα}c_{kα}⁺c_{kα}$
where:
*   $ε_{kα}$ is the energy of an electron in lead α with momentum k.
*   $c_{kα}⁺, c_{kα}$ are the creation and annihilation operators for an electron in lead α with momentum k.

###   Tunneling Hamiltonian $H_{tunneling}$ :
    
$H_{tunneling} = \sum_{kα,i} (t_{kα,i} c_{kα}⁺d_i + t_{kα,i}^* d_i⁺c_{kα})$

where:
*   $c_{kα}⁺$ ($c_{kα}$) creates (annihilates) an electron in lead α with momentum k
*   $d_i⁺$ ($d_i$) creates (annihilates) an electron in dot state i
*   $t_{kα,i}$ is the tunneling amplitude

The tunneling amplitudes are transformed from single-particle to many-body basis using `rotate_Tba()` in `qmeq/leadstun.py`.

## Step 2: Calculate Transition Rates $W_{ji}$

The transition rates $W_{ji}$ are calculated using Fermi's Golden Rule:

$ W_{ji} = (2π/ħ) |⟨i|H_{tunneling}|j⟩|^2 ρ(E_i) f(E_i - E_j - μ) $

where:
*   $|i⟩$ and $|j⟩$ are the final and initial many-body states
*   $ρ(E_i)$ is absorbed into the tunneling amplitudes (wide-band limit)
*   $f(x)$ is the Fermi-Dirac distribution: $f(x) = 1/(exp(x) + 1)$
*   $E_i - E_j$ is the energy difference between final and initial states
*   $μ$ is the chemical potential of the lead

### Implementation Details:

The rates are calculated in several steps:

1. **Matrix Elements** (in `generate_fct`):
```python
# Calculate |⟨i|H_tunneling|j⟩|^2 in many-body basis
xcb = (Tba[l, b, c]*Tba[l, c, b]).real
```

2. **Energy Conservation** (in `func_pauli`):
```python
# Normalized energy for numerical stability
alpha = (Ecb-mu)/T  # Ecb = E_c-E_b
# Fermi factors for addition/removal
cur0 = fermi_func(alpha)      # f((E_c-E_b-μ)/T)
cur1 = 1-cur0                 # 1-f((E_c-E_b-μ)/T)
```

3. **Assembly of Rates** (in `generate_coupling_terms`):
```python
# Combine matrix elements and Fermi factors
fctm -= paulifct[l, ba, 1]  # Electron removal
fctp += paulifct[l, ba, 0]  # Electron addition
```

### Important Notes:

1. **Wide-band Limit:** The density of states ρ(E) is absorbed into tunneling amplitudes
2. **Energy Conservation:** Energy differences are properly accounted for in many-body basis
3. **Fermi Statistics:** 
   - Addition (N_i = N_j + 1): Uses f(E_i-E_j-μ)
   - Removal (N_i = N_j - 1): Uses 1-f(E_i-E_j-μ)

## Step 3: Construct the PME Matrix

The Pauli master equation describes the time evolution of the probabilities P_i:

$\frac{dP_i}{dt} = \sum_j (W_{ij}P_j - W_{ji}P_i)$

The kernel matrix L is constructed in `generate_kern()` where:
- Diagonal terms: $L_{ii} = -\sum_j W_{ji}$ (total out-rate)
- Off-diagonal terms: $L_{ij} = W_{ij}$ (transition rate in)

This ensures probability conservation: $\sum_i \frac{dP_i}{dt} = 0$

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
