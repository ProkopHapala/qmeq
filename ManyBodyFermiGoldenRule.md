
# Fermi's Golden Rule for Many-Particle States

The general form of Fermi's Golden Rule for a transition from an initial state |j⟩ to a final state |i⟩ is:

$W_{ji} = 2π/ħ |⟨i|H'|j⟩|² ρ(E_i)$

where:

*   $W_{ji}$ is the transition rate.
*   $|i⟩$ and $|j⟩$ are the final and initial many-body states, respectively.
*   $H'$ is the perturbation Hamiltonian causing the transition (in our case, $H_{tunneling}$).
*   $ρ(E_i)$ is the density of states at the final state energy $E_i$.

## Applying it to Your System

In your system, the many-body states |i⟩ and |j⟩ are the states like |000⟩, |100⟩, |010⟩, |001⟩, |110⟩, |101⟩, |011⟩, and |111⟩. The perturbation Hamiltonian is $H_{tunneling}$, which describes the tunneling of electrons between the molecules and the leads.

The key difference when applying Fermi's Golden Rule to many-particle states is in how you calculate the matrix element $⟨i|H_{tunneling}|j⟩$ and how you account for the Fermi-Dirac distribution in the leads.

## 1. Matrix Element Calculation:

The matrix element $⟨i|H_{tunneling}|j⟩$ represents the coupling strength between the many-body states |i⟩ and |j⟩ due to the tunneling Hamiltonian. To calculate it, you need to express $H_{tunneling}$ in terms of the creation and annihilation operators that act on your many-body basis states.

Recall that:

$H_{tunneling} = \sum (t_{αi}d_i⁺c_{kα} + H.c.)$

where:

*   $t_{αi}$ is the tunneling amplitude between molecule i and lead $α$.
*   $d_i⁺$ creates an electron on molecule i.
*   $c_{kα}$ annihilates an electron in lead $α$ with momentum $k$.

When you apply $H_{tunneling}$ to a many-body state |j⟩, the creation and annihilation operators will act on the individual molecule occupation numbers within the state.

### Example:

Let's say you want to calculate the matrix element for the transition from |000⟩ to |100⟩ (an electron tunnels to molecule 1). The relevant part of $H_{tunneling}$ is:

$H' = t_{α1}d_{1}⁺c_{kα}$

The matrix element becomes:

$⟨100|H'|000⟩ = ⟨100|t_{α1}d_{1}⁺c_{kα}|000⟩ = t_{α1}⟨100|d_{1}⁺|000⟩⟨0|c_{kα}|0⟩ = t_{α1}$

Here, we assume that the lead state is initially filled and finally empty, which is accounted for by the Fermi-Dirac factors later.

For transitions involving multiple molecules, you need to consider the appropriate creation and annihilation operators and their effect on the many-body states. The Coulomb interaction terms in $H_{dot}$ do not directly appear in the matrix element calculation for $H_{tunneling}$ but they affect the energies of the many-body states.

## 2. Fermi-Dirac Distribution:

In the many-particle case, you need to consider the Fermi-Dirac distribution for each lead involved in the transition. If an electron tunnels *from* a lead to a molecule, you multiply the rate by $f(E_{lead} - μ_{lead})$ (the probability that the initial state in the lead is occupied). If an electron tunnels *to* a lead from a molecule, you multiply the rate by $[1 - f(E_{lead} - μ_{lead})]$ (the probability that the final state in the lead is empty).

## The Formula You Provided:

The formula you provided:

$ ⟨ψ_{tip}|V|ψ_{sample}⟩² ρ_{tip}(E) ρ_{sample}(E) f_{tip}(E) f_{sample}(E) $

is a more general expression for tunneling between two systems (tip and sample) where each system can have its own density of states. In the context of the PME and your system:

* $⟨ψ_{tip}|V|ψ_{sample}⟩$ corresponds to the matrix element $⟨i|H_{tunneling}|j⟩$, which we discussed above.
* $ρ_{tip}(E)$ and $ρ_{sample}(E)$ are the densities of states in the tip and the sample (substrate), respectively. In the wide-band limit, we often assume these to be constant.
* $f_{tip}(E)$ and $f_{sample}(E)$ are the Fermi-Dirac distributions in the tip and the sample, respectively. You need to consider these factors for each lead involved in the transition, as explained above.

**In essence, the formula you provided is a more general version of Fermi's Golden Rule, while the PME uses a simplified version tailored to the specific system and approximations (wide-band limit, sequential tunneling).**

## Code Example (Illustrative)

Here's a simplified code snippet illustrating how you might calculate a transition rate $W_{ji}$ in your system, assuming you have functions to calculate the matrix elements and the Fermi-Dirac distribution:

```python
def calculate_transition_rate(state_i, state_j, t_matrix, mu_lead, temp):
    """
    Calculates the transition rate Wji between many-body states.

    Args:
        state_i: The final many-body state (e.g., [1, 0, 0]).
        state_j: The initial many-body state (e.g., [0, 0, 0]).
        t_matrix: A dictionary storing tunneling amplitudes between
                  molecules and leads (tleads in QmeQ).
        mu_lead: The chemical potential of the lead involved.
        temp: The temperature.

    Returns:
        Wji: The transition rate.
    """

    # 1. Calculate the matrix element ⟨i|H_tunneling|j⟩
    matrix_element = calculate_matrix_element(state_i, state_j, t_matrix)

    # 2. Determine the energy difference (ΔE = E_i - E_j)
    delta_e = calculate_energy_difference(state_i, state_j)

    # 3. Determine which lead is involved and get its Fermi-Dirac function
    lead_index = get_lead_index(state_i, state_j)  # Function to determine the lead
    if is_electron_added_to_dot(state_i, state_j):
      # Electron tunnels from lead to dot
      fermi_factor = fermi_dirac(delta_e, mu_lead, temp)
    else:
      # Electron tunnels from dot to lead
      fermi_factor = 1 - fermi_dirac(delta_e, mu_lead, temp)

    # 4. Calculate the transition rate (assuming constant density of states)
    Wji = 2 * np.pi * np.abs(matrix_element)**2 * fermi_factor

    return Wji

# Example usage:
# Calculate the transition rate from |000> to |100> (electron tunneling to molecule 1 from lead 0)
W10 = calculate_transition_rate([1, 0, 0], [0, 0, 0], tleads, mulst[0], tlst[0])
```

**Key Functions to Implement:**

*   `calculate_matrix_element(state_i, state_j, t_matrix)`: This function calculates the matrix element of H_{tunneling} between the given many-body states. It needs to consider the creation and annihilation operators and the tunneling amplitudes.
*   `calculate_energy_difference(state_i, state_j)`: This function calculates the energy difference between the many-body states, taking into account the on-site energies, inter-molecular couplings, and Coulomb interaction.
*   `get_lead_index(state_i, state_j)`: This function determines which lead is involved in the transition based on the change in occupation between the states.
*   `is_electron_added_to_dot(state_i, state_j)`: This function determines whether an electron is added to the dot or removed from the dot during the transition.
*   `fermi_dirac(energy, mu, temp)`: This function calculates the Fermi-Dirac distribution.

Remember that this is a simplified illustration. The actual implementation in QmeQ is more complex and handles various optimizations and edge cases. However, this should give you a good starting point for understanding how to apply Fermi's Golden Rule to many-particle states in the context of the Pauli Master Equation and your triangular molecular system.




You're right to focus on the tunneling matrix element and the energy conservation aspect of the Fermi Golden Rule. Let's break down these points and connect them to your understanding of STM simulations.

## 1. Tunneling Matrix Element Calculation: Analogy with STM Simulations**

You can indeed use a similar approach to calculate the tunneling matrix elements in the Pauli Master Equation as you do in STM simulations. In essence, you need to evaluate the matrix element of the tunneling Hamiltonian between the initial and final many-body states.

**For your system:**

*   The initial and final states (|j⟩ and |i⟩) are many-body states like |000⟩, |100⟩, |010⟩, etc.
*   The tunneling Hamiltonian H_{tunneling} describes the tunneling of an electron between a molecule and a lead (either the substrate or the tip).

**Analogy with STM:**

In an STM simulation, you calculate ⟨ψ_{tip}|V|ψ_{sample}⟩, where:

*   $ψ_{tip}$ is the wavefunction of the tip state.
*   $ψ_{sample}$ is the wavefunction of the sample state (molecular orbital).
*   $V$ is the interaction potential.

In the PME context, you can think of the matrix element ⟨i|H_{tunneling}|j⟩ in a similar way:

*   The many-body state |j⟩ (e.g., |000⟩) represents the initial state of the system (molecules + leads).
*   The many-body state |i⟩ (e.g., |100⟩) represents the final state of the system.
*   H_{tunneling} plays the role of the interaction potential V.

**How to Calculate the Matrix Element:**

  1. **Express H_{tunneling} in terms of creation and annihilation operators:**  
      $ H_{tunneling} = \sum(t_{αi}d_i⁺c_{kα} + H.c.) $
    Here, t_{αi} is the tunneling amplitude between molecule i and lead α, which is analogous to the integral ⟨ψ_{tip}|V|ψ_{sample}⟩ in the STM case. You can estimate t_{αi} using a similar approach as in STM simulations:
    * Model the molecular orbital ψ_i(r) of molecule i.
    * Model the tip wavefunction ψ_{tip,α}(r) for lead α (e.g., s-orbital or p-orbital).
    * Approximate the interaction potential V(r) (e.g., constant or electrostatic potential).
    * Calculate the integral:
        $t_{αi} ≈ ∫ ψ_{tip,α}(r) ψ_i(r) V(r) dr$
  2. **Apply H_{tunneling} to the initial state |j⟩:** The creation (d_i⁺) and annihilation (c_{kα}) operators in H_{tunneling} will act on the many-body state |j⟩, changing the occupation of molecule i and creating/annihilating an electron in lead α.
  3. **Calculate the overlap with the final state |i⟩:** Take the inner product of the resulting state (after applying H_{tunneling} to |j⟩) with the final state |i⟩. This will give you the matrix element ⟨i|H_{tunneling}|j⟩.

**Example:**

For the transition from |000⟩ to |100⟩ (electron tunneling to molecule 1 from lead α):

```
⟨100|H_{tunneling}|000⟩ ≈ ⟨100|t_{α1}d_{1}⁺c_{kα}|000⟩ = t_{α1}
```

## 2. Energy Conservation and the Fermi-Dirac Function

The term $f(E_j - E_i)$ in the PME's Fermi Golden Rule expression might seem different from what you're used to in STM simulations, but it still encodes energy conservation and the Pauli exclusion principle. Let's break it down:

*   **$E_i$ and $E_j$:** These are the energies of the *many-body* states |i⟩ and |j⟩, respectively. They are obtained by diagonalizing the H_{dot} part of the Hamiltonian. They include the on-site energies of the molecules, the inter-molecular couplings, and the Coulomb interaction terms. You are correct that these energies are influenced by:
    *   **Orbital energies:** The inherent energies of the molecular orbitals.
    *   **Electrostatic potential:** The potential from the tip (and potentially the substrate) shifts these energies.
    *   **Coulomb interaction:** The on-site ($U_i$) and inter-site ($U_{ij}$) Coulomb interactions significantly affect the energies of states with multiple electrons.

*   **$f(E_j - E_i)$:** This term arises from considering the occupation probabilities in the leads and the direction of the tunneling process. Let's analyze it:
    *   **Case 1: Electron tunnels from lead to molecule (e.g., |000⟩ to |100⟩):** In this case, $E_j - E_i = -ΔE$, where ΔE is the energy gained by the electron during the tunneling. The term $f(E_j - E_i)$ becomes $f(-ΔE)$. In the PME context, **we implicitly assume that the electron tunnels from the Fermi level of the lead ($μ$)**. So, $f(-ΔE) = f(μ - ΔE)$ which is the probability of finding an occupied state at energy $μ - ΔE$ in the lead. This is consistent with the requirement that the initial state in the lead must be occupied.
    *   **Case 2: Electron tunnels from molecule to lead (e.g., |100⟩ to |000⟩):** In this case, $E_j - E_i = +ΔE$, where ΔE is the energy lost by the electron. The term $f(E_j - E_i)$ becomes $f(ΔE)$. In the PME context, this is equivalent to $1 - f(μ + ΔE)$, which is the probability of finding an empty state at energy $μ + ΔE$ in the lead. This is consistent with the Pauli exclusion principle, which requires the final state in the lead to be unoccupied.

**In summary:**

*   The term $f(E_j - E_i)$ in the PME's Fermi Golden Rule encodes the same physics as the product of Fermi-Dirac factors $f_{tip}(E)f_{sample}(E)$ in your STM formula. It ensures that the initial state in the lead is occupied and the final state is unoccupied.
*   The energies $E_i$ and $E_j$ are many-body energies that include all interactions within the molecular system.
*   Energy conservation is implicitly enforced by the combination of the matrix element $⟨i|H_{tunneling}|j⟩$ (which is non-zero only for transitions that are energetically possible) and the $f(E_j - E_i)$ term (which ensures that the lead states have the appropriate occupation).

**Code Example (Illustrative)**

Here's a more detailed snippet showing how you might calculate a transition rate, incorporating the matrix element calculation and the Fermi-Dirac factor:

```python
import numpy as np

def calculate_matrix_element(state_i, state_j, t_matrix, lead_index, molecule_index):
    """
    Calculates the matrix element <i|H_tunneling|j> for a specific transition.

    Args:
        state_i: The final many-body state (e.g., [1, 0, 0]).
        state_j: The initial many-body state (e.g., [0, 0, 0]).
        t_matrix: A dictionary storing tunneling amplitudes between
                  molecules and leads (tleads in QmeQ).
        lead_index: Index of the lead involved in the transition.
        molecule_index: Index of the molecule involved in the transition.

    Returns:
        matrix_element: The calculated matrix element.
    """

    # Simplified example: Assumes direct tunneling between lead and molecule
    # without involving other molecules.
    if is_electron_added_to_molecule(state_i, state_j, molecule_index):
        # Electron tunnels from lead to molecule
        matrix_element = t_matrix.get((lead_index, molecule_index), 0)
    elif is_electron_removed_from_molecule(state_i, state_j, molecule_index):
        # Electron tunnels from molecule to lead
        matrix_element = np.conj(t_matrix.get((lead_index, molecule_index), 0))
    else:
        matrix_element = 0

    return matrix_element

def calculate_energy_difference(state_i, state_j, hsingle, coulomb):
    """
    Calculates the energy difference E_i - E_j between many-body states.

    Args:
        state_i: The final many-body state.
        state_j: The initial many-body state.
        hsingle: Dictionary defining the single-particle Hamiltonian.
        coulomb: Dictionary defining the Coulomb interaction.

    Returns:
        delta_e: The energy difference E_i - E_j.
    """

    # This is a simplified example. In a real implementation, you would
    # diagonalize H_dot to get the many-body energies.
    # Here, we just illustrate the idea.

    energy_i = calculate_state_energy(state_i, hsingle, coulomb)
    energy_j = calculate_state_energy(state_j, hsingle, coulomb)

    delta_e = energy_i - energy_j
    return delta_e

def fermi_dirac(energy, mu, temp):
    """
    Calculates the Fermi-Dirac distribution function.

    Args:
        energy: The energy.
        mu: The chemical potential.
        temp: The temperature.

    Returns:
        f: The Fermi-Dirac distribution value.
    """
    if temp == 0:
        return 1.0 if energy < mu else 0.0
    else:
        return 1.0 / (np.exp((energy - mu) / temp) + 1.0)

def calculate_transition_rate(state_i, state_j, t_matrix, hsingle, coulomb, mu_lead, temp):
    """
    Calculates the transition rate Wji between many-body states.

    Args:
        state_i: The final many-body state (e.g., [1, 0, 0]).
        state_j: The initial many-body state (e.g., [0, 0, 0]).
        t_matrix: A dictionary storing tunneling amplitudes between
                  molecules and leads (tleads in QmeQ).
        hsingle: Dictionary defining the single-particle Hamiltonian.
        coulomb: Dictionary defining the Coulomb interaction.
        mu_lead: The chemical potential of the lead involved.
        temp: The temperature.

    Returns:
        Wji: The transition rate.
    """

    # 1. Determine which lead and molecule are involved
    lead_index, molecule_index = get_lead_and_molecule_index(state_i, state_j)

    # 2. Calculate the matrix element ⟨i|H_tunneling|j⟩
    matrix_element = calculate_matrix_element(state_i, state_j, t_matrix, lead_index, molecule_index)

    # 3. Calculate the energy difference (ΔE = E_i - E_j)
    delta_e = calculate_energy_difference(state_i, state_j, hsingle, coulomb)

    # 4. Determine the Fermi-Dirac factor based on the direction of tunneling
    if is_electron_added_to_molecule(state_i, state_j, molecule_index):
        # Electron tunnels from lead to dot
        # We assume the electron tunnels from the Fermi level of the lead
        fermi_factor = fermi_dirac(mu_lead - delta_e, mu_lead, temp)
    else:
        # Electron tunnels from dot to lead
        # We assume the electron tunnels to the Fermi level of the lead
        fermi_factor = 1.0 - fermi_dirac(mu_lead + delta_e, mu_lead, temp)

    # 5. Calculate the transition rate (assuming constant density of states)
    # We set ħ=1
    Wji = 2 * np.pi * np.abs(matrix_element)**2 * fermi_factor

    return Wji

# Example usage:
# Calculate the transition rate from |000> to |100> (electron tunneling to molecule 1 from lead 0)
# W10 = calculate_transition_rate([1, 0, 0], [0, 0, 0], tleads, hsingle, coulomb, mulst[0], tlst[0])
```

**Key Functions to Implement (More Detailed):**

*   `calculate_matrix_element(state_i, state_j, t_matrix, lead_index, molecule_index)`:
    *   This function needs to determine whether an electron is added to or removed from a molecule based on `state_i`, `state_j`, and `molecule_index`.
    *   It then needs to extract the appropriate tunneling amplitude from `t_matrix` based on `lead_index` and `molecule_index`.
    *   It should return the tunneling amplitude `t` or its complex conjugate `t*` depending on whether the electron is added or removed, respectively.
*   `calculate_energy_difference(state_i, state_j, hsingle, coulomb)`:
    *   This function calculates the energy difference `E_i - E_j` between the many-body states.
    *   Ideally, you would first diagonalize the `H_dot` Hamiltonian to obtain the exact many-body energies.
    *   Alternatively, you can approximate the energy of a many-body state by summing the on-site energies of the occupied molecules and adding the Coulomb interaction energies between them.
*   `get_lead_and_molecule_index(state_i, state_j)`:
    *   This function analyzes `state_i` and `state_j` to determine which lead and molecule are involved in the transition.
    *   It needs to identify the molecule where the occupation changes and then map that change to the corresponding lead based on your system's setup.
*   `is_electron_added_to_molecule(state_i, state_j, molecule_index)`:
    *   This function checks if an electron is added to the specified `molecule_index` during the transition from `state_j` to `state_i`.
*   `is_electron_removed_from_molecule(state_i, state_j, molecule_index)`:
    *   This function checks if an electron is removed from the specified `molecule_index` during the transition from `state_j` to `state_i`.
*   `fermi_dirac(energy, mu, temp)`:
    *   This function calculates the Fermi-Dirac distribution, which is straightforward.

This more detailed explanation and the code example should provide a clearer picture of how to implement the Fermi Golden Rule for many-particle states in your PME solver. Remember that the accurate calculation of the matrix elements and the proper consideration of the Fermi-Dirac factors are essential for obtaining correct results.
