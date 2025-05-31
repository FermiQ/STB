# Documentation for `src/hamiltionian.f90`

## Overview

The `Class_hamiltionian` module defines the electronic Hamiltonian for a crystal structure. It encapsulates the parameters and methods required to construct the Hamiltonian matrix for given k-points, calculate its eigenvalues and eigenvectors, and compute derivatives of the Hamiltonian with respect to k (which are essential for velocity operators and Berry curvature calculations). The module supports various types of Hamiltonians, including s-orbital and p-orbital models, nearest and second-nearest neighbor hopping, spin-orbit coupling (Rashba and intrinsic), local exchange splitting, and specialized models like the Haldane model and Kane-Mele model. It also handles MPI for parallel execution and reading configurations.

## Key Components

### `module Class_hamiltonian`

This is the main module.

### `type hamil`

A derived data type that holds all the parameters and procedures related to the Hamiltonian.

**Key Members (Parameters):**

-   `E_fermi (real(8), allocatable :: E_fermi(:))`: Array of Fermi levels.
-   `temp (real(8))`: Temperature used in the Fermi-Dirac distribution.
-   `E_s (real(8))`: Onsite energy for s-orbitals.
-   `E_A, E_B (real(8))`: Onsite energies for A and B sublattice sites (e.g., in honeycomb).
-   `E_p(3) (real(8))`: Onsite energies for p_x, p_y, p_z orbitals.
-   `Vss_sig (real(8))`: Nearest neighbor hopping for s-orbitals (sigma bond).
-   `Vpp_pi, Vpp_sig (real(8))`: Nearest neighbor hopping for p-orbitals (pi and sigma bonds).
-   `V2pp_pi, V2pp_sig (real(8))`: Second nearest neighbor hopping for p-orbitals.
-   `eta_soc (real(8))`: Parameter for p-state intrinsic spin-orbit coupling.
-   `t_2 (real(8))`: Amplitude for 2nd nearest neighbor hopping (Haldane model).
-   `phi_2 (real(8))`: Phase for 2nd nearest neighbor hopping (Haldane model).
-   `t_so (real(8))`: Rashba spin-orbit coupling strength.
-   `lambda (real(8))`: Local exchange splitting strength.
-   `HB1, HB2, HB_eta (real(8))`: Parameters for "Hongbin's model".
-   `lambda_KM (real(8))`: Parameter for the Kane-Mele term.
-   `gamma (real(8))`: Broadening factor (e.g., for Green's functions, conductivity).
-   `drop_Vx_layers(:), drop_Vy_layers(:) (real(8), allocatable)`: Specifies layers where hopping derivatives (velocities) in x or y are zeroed out.
-   `del_H (complex(8), allocatable :: del_H(:,:))`: Matrix for the derivative of H w.r.t. k.
-   `prefix (character(len=300))`: Output file prefix.
-   `nProcs (integer)`: Number of MPI processes.
-   `me (integer)`: MPI rank of the current process.
-   `num_orb (integer)`: Number of orbitals per atom (e.g., 1 for s, 3 for p).
-   `num_up (integer)`: Number of spin-up states (num_orb * num_atoms_in_UC).
-   `test_run (logical)`: Flag to perform unit tests.
-   `UC (type(unit_cell))`: Unit cell object.
-   `units (type(units))`: Units object for physical quantities.

**Key Procedures (Methods):**

-   **`init_hamil(cfg) result(self)`**: Constructor. Initializes the `hamil` object from a configuration (`cfg`). Reads parameters, broadcasts them using MPI.
-   **`Bcast_hamil(self)`**: Broadcasts Hamiltonian parameters to all MPI ranks.
-   **`setup_H(self, k, H)`**: Constructs the Hamiltonian matrix `H` for a given k-point `k`. This is a central routine that calls various `set_*` subroutines based on the enabled parameters.
-   **`calc_eigenvalues(self, k_list, eig_val)`**: Calculates eigenvalues `eig_val` for a list of k-points `k_list`. Uses LAPACK's `zheevd`.
-   **`calc_single_eigenvalue(self, k, eig_val)`**: Calculates eigenvalues for a single k-point.
-   **`calc_eig_and_velo(self, k, eig_val, del_kx, del_ky, pert_log)`**: Calculates eigenvalues and velocity matrices (derivatives of H) `del_kx`, `del_ky`. Can include perturbative corrections (`pert_log`).
-   **`set_fermi(self, cfg)`**: Sets the Fermi energy levels based on the configuration.
-   **`fermi_distr(self, E, n_ferm) result(ferm)`**: Calculates the Fermi-Dirac distribution value.
-   **`set_EigenE(self, H)`**: Adds onsite eigenenergies (E_s, E_A, E_B) to the Hamiltonian `H`.
-   **`set_p_energy(self, H)`**: Adds onsite p-orbital energies (E_p) to `H`.
-   **`set_hopping(self, k, H)`**: Adds nearest-neighbor hopping terms to `H`. Calls `set_hopp_mtx`.
-   **`set_snd_hopping(self, k, H)`**: Adds second-nearest-neighbor hopping terms to `H`. Calls `set_snd_hopp_mtx`.
-   **`set_haldane_hopping(self, k, H)`**: Adds Haldane model specific second-neighbor hopping terms.
-   **`set_rashba_SO(self, k, H)`**: Adds Rashba spin-orbit coupling terms.
-   **`set_SOC(self, H)`**: Adds intrinsic p-orbital spin-orbit coupling. Calls `set_small_SOC`.
-   **`set_loc_exch(self, H)`**: Adds local exchange splitting (Stoner magnetism). Calls `setup_Stoner_mtx`.
-   **`set_KaneMele_exch(self, k, H)`**: Adds Kane-Mele model terms.
-   **`set_hongbin_hopp(self, k, H)`**: Adds hopping terms for "Hongbin's model".
-   **`set_hongbin_SOC(self, H)`**: Adds SOC terms for "Hongbin's model".
-   **`set_hopp_mtx(self, R, hopp_mtx)`**: Sets up the orbital-space hopping matrix for a given connection vector `R` (for s or p orbitals).
-   **`set_snd_hopp_mtx(self, R, hopp_mtx)`**: Similar to `set_hopp_mtx` but for second-nearest neighbors (p-orbitals).
-   **`set_p_hopp_mtx(R, Vpp_sig, Vpp_pi, hopp_mtx)`**: Helper for p-orbital hopping matrix construction based on Slater-Koster parameters.
-   **`set_derivative_k(self, k, k_idx)`**: Calculates `self%del_H`, the derivative of the Hamiltonian with respect to `k(k_idx)`. Calls various `set_derivative_*` subroutines.
-   **`set_deriv_FD(self, k, k_idx, del_H)`**: Calculates `del_H` using finite differences (for testing or for terms without analytical derivatives).
-   **`set_derivative_hopping(self, k, k_idx)`**: Analytical derivative for nearest-neighbor hopping.
-   **`set_derivative_snd_hopping(self, k, k_idx)`**: Analytical derivative for second-nearest-neighbor hopping.
-   **`set_derivative_haldane_hopping(self, k, k_idx)`**: Analytical derivative for Haldane hopping.
-   **`set_derivative_rashba_so(self, k, k_idx)`**: Analytical derivative for Rashba SOC.
-   **`set_derivative_KM(self, k)`**: Analytical derivative for Kane-Mele term.
-   **`compare_derivative(self, k)`**: Compares analytical derivatives with finite difference calculations for testing.
-   **`calc_velo_mtx(self, k, derive_idx, eig_vec_mtx, ret)`**: Computes velocity operator matrix elements in the eigenbasis: `V = U* (dH/dk) U`, where U are eigenvectors.
-   **`calc_left_pert_velo_mtx`, `calc_right_pert_velo_mtx`**: Calculate velocity matrices with first-order perturbative corrections due to exchange splitting.
-   **`calc_exch_firstord`**: Calculates the first-order correction to the Hamiltonian due to changes in exchange field.
-   **`calc_berry_z(self, z_comp, eig_val, x_mtx, y_mtx)`**: Calculates the z-component of Berry curvature for each band.
-   **`calc_berry_diag`, `calc_berry_diag_sea`, `calc_berry_diag_surf`**: Calculate diagonal components of conductivity-like tensors (related to Berry curvature) using different formulations (e.g., Kubo-Streda, considering Fermi sea/surface terms).
-   **`calc_fac_diag`, `calc_fac_sea`, `calc_fac_surf`**: Helper functions to calculate energy-dependent factors in Berry curvature expressions.
-   **`drop_layer_derivative(self, k_idx)`**: Zeros out elements of `self%del_H` corresponding to specified layers.
-   **`z_layer_states(self) result(z)`**: Returns the z-coordinates of atomic layers.
-   **`free_ham(self)`**: Deallocates memory associated with the `hamil` object.
-   **`test_herm(H, tag, verbose)`**: Tests if a matrix `H` is Hermitian.

## Important Variables/Constants (within `type hamil`)

-   **`E_fermi(:)`**: Array storing Fermi energy values. Crucial for determining occupations and transport properties.
-   **`Vss_sig, Vpp_sig, Vpp_pi, V2pp_sig, V2pp_pi`**: Hopping parameters defining the band structure.
-   **`t_so`**: Rashba spin-orbit coupling strength.
-   **`eta_soc`**: Intrinsic p-orbital SOC strength.
-   **`lambda`**: Local exchange field strength (for magnetism).
-   **`lambda_KM`**: Kane-Mele SOC strength.
-   **`t_2, phi_2`**: Haldane model parameters.
-   **`gamma`**: Broadening parameter used in Green's functions and conductivity calculations.
-   **`num_orb`**: Number of orbitals per atom site.
-   **`num_up`**: Total number of spin-up basis states in the unit cell. The Hamiltonian matrix size is `(2 * num_up) x (2 * num_up)`.
-   **`del_H(:,:)`**: The derivative of the Hamiltonian, dH/dk, used to calculate velocities.
-   **`me, root`**: MPI rank and root process ID.
-   Constants like `sigma_x`, `sigma_y`, `sigma_z` (Pauli matrices), `i_unit` (imaginary unit) are likely used from the `Constants` module.

## Dependencies/Interactions

-   **`m_config`**: Used for reading configuration parameters (`CFG_t`, `CFG_get`, `CFG_get_size`).
-   **`output`**: Used for error messages and logging (`error_msg`).
-   **`Class_unit_cell`**: Provides the `unit_cell` type (`UC`) which describes the crystal lattice, atomic positions, and connectivity (neighbors). This is fundamental for constructing the Hamiltonian.
-   **`m_npy`**: Used for saving data in NumPy `.npy` format (e.g., `save_npy` for Fermi levels).
-   **`mpi`**: Standard MPI module for parallel communication (`MPI_COMM_WORLD`, `MPI_Bcast`, `MPI_Comm_rank`, `MPI_Comm_size`, etc.).
-   **`MYPI`**: Custom MPI utilities, possibly defining `MYPI_INT` or wrapper functions.
-   **`Constants`**: Provides physical constants (e.g., `boltzmann_const`, `PI`, Pauli matrices `sigma_x, sigma_y, sigma_z`, imaginary unit `i_unit`, `c_0`, `c_1`, `c_i`) and possibly numerical precision parameters (`eta_sq`, `pos_eps`).
-   **LAPACK/BLAS**: External libraries used for numerical linear algebra, specifically `zheevd` (or `zheev`) for diagonalizing Hermitian matrices to find eigenvalues and eigenvectors, and `zgemm` for matrix-matrix multiplication. These are implicitly used via Fortran's intrinsic interfaces or direct calls.

The module interacts heavily with the `Class_unit_cell` to understand the geometry and connections between atoms, which dictates the non-zero elements of the Hamiltonian. It receives k-points as input and produces eigenvalues, eigenvectors, and velocity operators as output, which are then used by other parts of the larger program (e.g., `Class_k_space`) for calculating physical observables.
