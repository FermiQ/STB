# Documentation for `src/k_space.F90`

## Overview

The `Class_k_space` module is responsible for calculations that involve sampling the Brillouin zone (k-space). This includes calculating band structures along specified paths, computing the Density of States (DOS), and evaluating Berry curvature related quantities like the Anomalous Hall Conductivity (AHC) and orbital magnetization. It handles the generation of k-point grids, manages MPI parallelization for k-point loops, and interfaces with the `Class_hamiltonian` module to obtain eigenvalues and velocity operators at each k-point. The module also implements adaptive k-mesh refinement strategies for accurate integration of Berry curvature.

## Key Components

### `module Class_k_space`

This is the main module.

### `type k_space`

A derived data type that holds parameters and procedures related to k-space calculations.

**Key Members (Parameters):**

-   `new_k_pts (real(8), allocatable :: new_k_pts(:,:))`: Array storing newly generated k-points for a calculation step.
-   `all_k_pts (real(8), allocatable :: all_k_pts(:,:))`: Array storing all k-points accumulated during adaptive refinement.
-   `int_DOS (real(8), allocatable :: int_DOS(:))`: Integrated Density of States.
-   `E_DOS (real(8), allocatable :: E_DOS(:))`: Energy grid for DOS calculation.
-   `DOS_gamma (real(8))`: Broadening parameter ($\Gamma$) for DOS calculations.
-   `DOS_lower, DOS_upper (real(8))`: Lower and upper energy bounds for DOS calculation.
-   `temp (real(8))`: Temperature used in Fermi-Dirac distribution (primarily for Berry quantities).
-   `DOS_num_k_pts (integer)`: Number of k-points per dimension for DOS calculations on a regular grid.
-   `berry_num_k_pts (integer)`: Initial number of k-points per dimension for Berry curvature calculations.
-   `ACA_num_k_pts (integer)`: Number of k-points for Anomalous Hall Conductivity (ACA) calculations (likely a specific type of Berry calculation).
-   `num_DOS_pts (integer)`: Number of points on the energy grid for DOS.
-   `num_k_pts (integer)`: Number of k-points per segment for band structure paths.
-   `num_plot_pts (integer)`: Number of points per dimension for plotting Berry curvature.
-   `nProcs (integer)`: Number of MPI processes.
-   `me (integer)`: MPI rank of the current process.
-   `berry_iter (integer)`: Maximum number of iterations for adaptive k-mesh refinement.
-   `kpts_per_step (integer)`: Number of new k-points to add per refinement step per MPI process.
-   `k_shift(3) (real(8))`: Shift applied to the k-space grid.
-   `berry_conv_crit (real(8))`: Convergence criterion for Berry integration.
-   `weights (real(8), allocatable :: weights(:))`: Integration weights for each k-point (e.g., Voronoi cell area / number of triangle vertices).
-   `elem_nodes (integer, allocatable :: elem_nodes(:,:))`: Connectivity list for k-point triangulation (defines triangles).
-   `refine_weights(:), refine_weights_surf(:), refine_weights_sea(:) (real(8), allocatable)`: Weights used to guide adaptive mesh refinement, based on Berry curvature contributions.
-   `k1_param(:), k2_param(:) (real(8), allocatable)`: Parameters defining k-paths or grids (e.g., start/end points, number of points).
-   `filling (character(len=300))`: String indicating how k-points are generated (e.g., "path_rel", "path_abs", "grid").
-   `prefix (character(len=300))`: Output file prefix.
-   `chosen_weights (character(len=300))`: Which quantity to use for adaptive refinement weights ("hall" or "orbmag").
-   `ada_mode (character(len=6))`: Mode for adaptive refinement ("area", "weight", "mixed").
-   `berry_component (character(len=6))`: Specifies which component of Berry curvature or related tensor to calculate (e.g., "xy", "xx", "xysurf").
-   `perform_pad (logical)`: Whether to pad the k-point list to ensure an even distribution across MPI processes.
-   `calc_hall, calc_hall_diag, calc_orbmag (logical)`: Flags to enable calculation of Hall conductivity (standard Kubo), diagonal Hall conductivity (Kubo-Streda like), and orbital magnetization.
-   `test_run (logical)`: Flag to perform unit tests.
-   `berry_safe (logical)`: Flag to save intermediate Berry curvature data.
-   `ham (type(hamil))`: The Hamiltonian object.
-   `units (type(units))`: Units object.

**Key Procedures (Methods):**

-   **`init_k_space(cfg) result(self)`**: Constructor. Initializes `k_space` object from configuration `cfg`. Sets up parameters and broadcasts them.
-   **`Bcast_k_space(self)`**: Broadcasts `k_space` parameters to all MPI ranks.
-   **`calc_and_print_band(self)`**: Calculates and saves band structure.
    -   Calls `setup_k_path_rel`, `setup_k_path_abs`, or `setup_k_grid` based on `self%filling`.
    -   Distributes k-points among MPI processes.
    -   Calls `self%ham%calc_eigenvalues`.
    -   Gathers and saves k-points and eigenvalues.
-   **`calc_and_print_dos(self)`**: Calculates and saves DOS and projected DOS (PDOS).
    -   Sets up a k-grid using `setup_inte_grid_para` or `setup_inte_grid_hex`.
    -   Calls `calc_pdos`.
    -   Saves DOS, PDOS, and integrated DOS.
-   **`calc_pdos(self, E, PDOS)`**: Calculates projected density of states.
    -   Loops over k-points, calculates eigenvalues and eigenvectors using `zheevd`.
    -   Projects eigenvectors onto atomic orbitals and broadens with `lorentzian`.
-   **`calc_berry_quantities(self, pert_log)`**: Main routine for calculating Berry curvature related quantities (AHC, orbital magnetization).
    -   Calls `setup_berry_inte_grid` for initial k-mesh.
    -   Iteratively refines the k-mesh:
        -   `calc_new_berry_points`: Calculates eigenvalues, Berry curvature ($\Omega_z$), and orbital magnetization moments ($Q_L, Q_{IC}$) at new k-points using `self%ham%calc_eig_and_velo`.
        -   `append_kpts`, `append_eigval`, `append_quantity`: Adds new data to global lists.
        -   `integrate_hall`, `integrate_hall_surf`, `integrate_hall_sea`, `integrate_orbmag`: Performs integration over the current k-mesh using `set_weights_ksp` (based on triangulation via `run_triang`).
        -   `process_hall`, `process_hall_surf`, `process_orbmag`: Checks convergence and saves iteration data.
        -   `set_hall_weights`, `set_orbmag_weights`: Determines weights for next refinement step.
        -   `add_kpts_iter`: Selects new k-points based on `ada_mode` and weights.
    -   Calls `finalize_hall`, `finalize_hall_surf`, `finalize_orbmag` to save final results.
-   **`setup_k_path_rel(self)` / `setup_k_path_abs(self)`**: Generates k-points along a path defined by relative (to reciprocal lattice vectors) or absolute coordinates.
-   **`setup_k_grid(self)`**: Generates a regular 2D grid of k-points.
-   **`setup_inte_grid_para(self, n_k, padding)`**: Sets up a parallelogram-shaped k-grid based on reciprocal lattice vectors.
-   **`setup_inte_grid_hex(self, n_k)`**: Sets up a k-grid suitable for hexagonal Brillouin zones.
-   **`setup_berry_inte_grid(self)`**: Calls appropriate grid setup for Berry calculations.
-   **`integrate_hall(self, kidx_all, omega_z_all, eig_val_all, hall)`**: Integrates $\Omega_z \cdot f(E)$ to get Hall conductivity.
-   **`integrate_hall_surf(self, kidx_all, omega_z_all, hall)` / `integrate_hall_sea(...)`**: Integrates components for diagonal conductivity (Fermi surface and Fermi sea terms).
-   **`integrate_orbmag(self, Q_kidx_all, Q_L_all, Q_IC_all, orb_mag, orbmag_L, orbmag_IC)`**: Integrates local ($Q_L$) and itinerant circulation ($Q_{IC}$) moments to get orbital magnetization.
-   **`set_weights_ksp(self)`**: Calculates integration weights for k-points based on the area of surrounding triangles from `run_triang`.
-   **`lorentzian(self, x) result(lor)`**: Calculates Lorentzian broadening function.
-   **`find_fermi(self, cfg)`**: Determines Fermi energy based on target filling from integrated DOS.
-   **`calc_ACA(self)`**: Calculates Anomalous Hall Conductivity using a specific formulation (likely related to orbital moments of p-orbitals). Calls `calc_ACA_singleK`.
-   **`plot_omega_square(self)`**: Generates data for plotting Berry curvature over a 2D k-grid.
-   **`free_ksp(self)`**: Deallocates memory.

## Important Variables/Constants

-   **`DOS_gamma (real(8))`**: Lorentzian broadening width for DOS calculations.
-   **`berry_num_k_pts (integer)`**: The initial density of the k-mesh for Berry calculations. The mesh is adaptively refined from this starting point.
-   **`k_shift(3) (real(8))`**: A vector used to shift the entire k-mesh. Useful for sampling around specific points or avoiding high-symmetry points if necessary.
-   **`filling (character(len=300))`**: Controls how k-points are generated for band structure (path, grid).
-   **`prefix (character(len=300))`**: Used for naming output files.
-   **`elem_nodes (integer, allocatable :: elem_nodes(:,:))`**: Stores the triangulation of the k-mesh, crucial for integration.
-   **`weights (real(8), allocatable :: weights(:))`**: Stores the integration weight for each k-point.
-   **`ada_mode (character(len=6))`**: Determines the strategy for adaptive k-mesh refinement (e.g., 'area' focuses on largest triangles, 'weight' focuses on largest Berry curvature contributions).
-   `PI`, `deg_30`, `deg_60`, `speed_of_light`, `eta`, `pos_eps` (likely from `Constants` or `Class_helper`): Mathematical and physical constants.

## Dependencies/Interactions

-   **`m_config`**: For reading configuration parameters (`CFG_t`, `CFG_get`, `create_dir`).
-   **`m_npy`**: For saving data in NumPy `.npy` format (`save_npy`).
-   **`mpi`**: Standard MPI module for parallelization (`MPI_COMM_WORLD`, `MPI_Bcast`, `MPI_Reduce`, `MPI_Gatherv`, etc.).
-   **`Class_hamiltonian`**: Crucial dependency. The `k_space` type contains a `hamil` object (`self%ham`). This is used to:
    -   Get Hamiltonian parameters.
    -   Calculate eigenvalues (`self%ham%calc_eigenvalues`, `self%ham%calc_eig_and_velo`).
    -   Calculate velocity operators (dH/dk via `self%ham%calc_eig_and_velo`, which calls `self%ham%set_derivative_k` and `self%ham%calc_velo_mtx`).
    -   Calculate Berry curvature components (`self%ham%calc_berry_z`, `self%ham%calc_berry_diag_surf`, `self%ham%calc_berry_diag_sea`).
    -   Get Fermi-Dirac distribution (`self%ham%fermi_distr`).
    -   Get unit cell information (`self%ham%UC%...`).
-   **`Class_helper`**: Provides various helper functions:
    -   `my_section`, `sections`: For distributing loop iterations among MPI processes.
    -   `linspace`: For generating evenly spaced points.
    -   `cross_prod`, `my_norm2`: Vector operations.
    -   `qargsort`: For sorting and getting indices.
    -   `find_list_idx`: To find an element in a list.
    -   `append_kidx`, `append_eigval`, `append_quantity`: Utilities for concatenating arrays during adaptive refinement.
    -   Possibly constants like `deg_30`, `deg_60`, `PI`, `eta`, `pos_eps`.
-   **`MYPI`**: Custom MPI utilities, possibly for error checking (`check_ierr`) or custom types (`MYPI_INT`).
-   **`ieee_arithmetic`**: For checking for NaN values (`ieee_is_nan`).
-   **`run_triang` (interface)**: External library/subroutine for Delaunay triangulation of the k-points to define `elem_nodes`. This is fundamental for the integration scheme.
-   **LAPACK**: Implicitly used via `zheevd` (called within `calc_pdos` and `self%ham%calc_eig_and_velo`, etc.) for eigenvalue/eigenvector computations.

The module orchestrates complex calculations by combining k-space sampling strategies with Hamiltonian evaluations and subsequent post-processing (integration, summation) to derive physical quantities.
