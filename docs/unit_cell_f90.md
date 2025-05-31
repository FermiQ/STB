# Documentation for `src/unit_cell.f90`

## Overview

The `Class_unit_cell` module defines the structure and magnetic configuration of the unit cell in a crystalline material. It handles the geometry of the lattice, including lattice vectors and reciprocal lattice vectors, the positions and types of atoms within the unit cell, and the connectivity between these atoms (nearest neighbors, second nearest neighbors). A significant part of this module is dedicated to setting up various magnetic textures, such as ferromagnetic, antiferromagnetic, spiral, and skyrmion configurations.

## Key Components

### `module Class_unit_cell`

This is the main module.

### `enum, bind(c)` for Connection Types

Defines enumerators for different types of neighbor connections:
-   `nn_conn`: Nearest Neighbor connection.
-   `snd_nn_conn`: Second Nearest Neighbor connection.
-   `oop_conn`: Out-of-plane connection (though primarily 2D systems seem to be handled).

### `type unit_cell`

A derived data type that encapsulates all information about the unit cell.

**Key Members (Parameters):**

-   `lattice(2, 2) (real(8))`: Real-space lattice vectors (2 vectors in 2D).
-   `rez_lattice(2, 2) (real(8))`: Reciprocal lattice vectors.
-   `num_atoms (integer)`: Number of non-redundant atoms in the unit cell.
-   `atom_per_dim (integer)`: Number of atoms along a characteristic dimension (e.g., radius for hexagonal cells, side length for square cells).
-   `nProcs, me (integer)`: MPI process count and rank.
-   `n_wind (integer)`: Winding number for certain magnetic textures (e.g., skyrmions, linear rotation spirals).
-   `wavevector (integer, allocatable :: wavevector(:))`: Wavevector `q` for spiral magnetic orders, typically in units of reciprocal lattice vectors.
-   `lattice_constant (real(8))`: The fundamental lattice spacing.
-   `eps (real(8))`: Tolerance for positional accuracy when identifying atoms or neighbors.
-   `ferro_phi, ferro_theta (real(8))`: Spherical angles ($\phi, \theta$) for uniform ferromagnetic alignment.
-   `cone_angle (real(8))`: Cone angle for conical spiral magnetic structures.
-   `anticol_phi(:), anticol_theta(:) (real(8), allocatable)`: Arrays of angles for defining antiferromagnetic or more complex collinear/non-collinear structures, often used for defining base states for perturbations or spirals.
-   `axis_phi, axis_theta (real(8))`: Spherical angles defining the rotation axis for spiral magnetism.
-   `atan_factor, atan_pref, skyrm_middle, dblatan_dist, dblatan_pref (real(8))`: Parameters for defining skyrmion profiles using `atan` or double `atan` functions for the magnetization tilt.
-   `atoms (type(atom), dimension(:), allocatable)`: Array of `atom` objects that constitute the unit cell. Each `atom` object (from `Class_atom`) stores its position, magnetic moment, site type, and neighbor information.
-   `units (type(units))`: Object storing physical units.
-   `uc_type (character(len=25))`: String identifier for the unit cell type (e.g., "square_2d", "honey_2d", "honey_line", "file_square").
-   `mag_type (character(len=25))`: String identifier for the type of magnetic ordering.
-   `spiral_type (character(len=25))`: Specific type of spiral (e.g., "anticol", "cone").
-   `mag_file (character(len=300))`: Path to a file defining a custom magnetic configuration.
-   `molecule (logical)`: If true, treats the system as an isolated molecule (likely affects boundary conditions or k-space interpretation).
-   `test_run (logical)`: Flag to run internal consistency tests.
-   `pert_log (logical)`: Flag indicating if Berry curvature calculations should use first-order perturbation theory (related to magnetic perturbations).

**Key Procedures (Methods):**

-   **`init_unit(cfg) result(self)`**: Constructor. Initializes the `unit_cell` object from a configuration `cfg`.
    -   Reads general parameters (epsilon, magnetic types, lattice constant, etc.).
    -   Broadcasts parameters using `Bcast_UC`.
    -   Calls a specific `init_unit_TYPE` subroutine based on `uc_type`.
    -   Calculates reciprocal lattice vectors.
-   **`Bcast_UC(self)`**: Broadcasts unit cell parameters to all MPI ranks.
-   **`init_unit_square(self)`**: Initializes a square 2D unit cell. Calls `setup_square` and then `setup_gen_conn`. Sets magnetic order based on `mag_type`.
-   **`init_unit_honey_hexa(self)`**: Initializes a hexagonal unit cell with honeycomb lattice structure. Uses `make_hexagon` and `setup_honey`, then `setup_gen_conn` and `set_honey_snd_nearest`. Sets magnetic order.
-   **`init_unit_honey_line(self)`**: Initializes a line-shaped (ribbon-like) unit cell with honeycomb structure. Uses `make_honeycomb_line`, `setup_honey`, `setup_gen_conn`, and `set_honey_snd_nearest_line`. Sets magnetic order.
-   **`init_file_square(self)`**: Initializes a unit cell from a file, typically defining atomic positions and magnetic moments.
-   **`setup_square(self)`**: Populates atomic positions for a square lattice.
-   **`setup_honey(self, hexagon, site_type)`**: Populates atomic positions for a honeycomb lattice based on pre-calculated `hexagon` coordinates and `site_type` (A or B).
-   **`make_hexagon(self, hexagon, site_type)`**: Generates the atomic coordinates and site types for a hexagonal supercell of a honeycomb lattice.
-   **`make_honeycomb_line(self, line, site_type)`**: Generates atomic coordinates and site types for a honeycomb ribbon.
-   **`setup_gen_conn(self, conn_mtx, conn_type, transl_mtx)`**: General routine to find and store neighbor information for each atom.
    -   `conn_mtx`: Defines relative vectors to potential neighbors.
    -   `conn_type`: Specifies the type of connection (nn_conn, snd_nn_conn).
    -   `transl_mtx`: Defines translation vectors to adjacent unit cells for periodic boundary conditions.
    -   Calls `gen_find_neigh` for each potential connection.
-   **`gen_find_neigh(self, start, conn, transl_mtx) result(neigh)`**: Finds the index of a neighboring atom given a starting atom, a connection vector, and periodic translation vectors. Uses `in_cell`.
-   **`in_cell(self, start, conn) result(idx)`**: Checks if the atom at `start + conn` is within the current list of atoms and returns its index.
-   **`set_honey_snd_nearest(self)` / `set_honey_snd_nearest_line(self)`**: Specifically sets up second-nearest neighbor connections for honeycomb structures.
-   **Magnetic Configuration Methods:**
    -   `set_mag_ferro(self)`: Sets all atomic magnetic moments to a uniform ferromagnetic state.
    -   `set_mag_anticol(self)`: Sets moments based on `self%anticol_phi` and `self%anticol_theta`, allowing for various collinear/non-collinear antiferromagnetic or complex states. Can also set a base collinear state if `pert_log` is true.
    -   `set_mag_random(self)`: Sets random magnetic moments.
    -   `set_mag_x_spiral_square(self)`: Sets a 1D spiral along x for a square lattice.
    -   `set_mag_linrot_1D_spiral(self, center, UC_l)`: Sets a 1D linear spiral/spin-wave. `m0_A` and `m0_B` define the base magnetic state(s) being rotated. The rotation angle depends on the projection of the atom's position onto the `wavevector`.
    -   `set_mag_linrot_1D_spiral_honey(self)`: Wrapper for 1D spiral on honeycomb, choosing between `anticol` or `cone` base states for the spiral.
    -   `set_mag_linrot_1D_spiral_m0_anticol(self)` / `set_mag_linrot_1D_spiral_m0_cone(self)`: Initialize `m0_A` and `m0_B` (base moments for the spiral) for anticollinear or conical configurations respectively.
    -   `set_mag_linrot_skyrm(self, center, radius)`: Sets a skyrmion magnetic texture where moments rotate based on their position relative to `center` and `radius`.
    -   `set_mag_atan_skyrm(self, center, radius)` / `set_mag_dblatan_skyrm(self, center, radius)`: Skyrmion profiles defined by `atan` or double `atan` functions for smoother or more controlled domain walls.
    -   `set_mag_linrot_skrym_square(self)` / `set_mag_linrot_skrym_honey(self)`: Wrappers for skyrmions on square or honeycomb lattices.
-   `save_unit_cell(self, folder)`: Saves unit cell information (positions, magnetic moments, lattice vectors, neighbor lists) to files using `m_npy`.
-   `calc_area(self) result(area)`: Calculates the 2D area of the unit cell.
-   `free_uc(self)`: Deallocates memory associated with the unit cell, including its atoms.
-   `run_tests(self)`: Performs internal consistency checks, e.g., for MPI broadcasted values.

## Important Variables/Constants (within `type unit_cell`)

-   **`lattice(2,2)`**: The two 2D real-space lattice vectors defining the unit cell.
-   **`rez_lattice(2,2)`**: The two 2D reciprocal lattice vectors.
-   **`num_atoms`**: Total number of atoms in the defined unit cell.
-   **`uc_type (character)`**: String specifying the lattice type (e.g., "square_2d", "honey_2d"). This determines which initialization routines are called.
-   **`mag_type (character)`**: String specifying the magnetic configuration (e.g., "ferro", "1Dspiral", "lin_skyrm"). This determines which `set_mag_TYPE` routine is called.
-   **`atoms(:) (type(atom))`**: The array holding all atomic information. The properties of these `atom` objects (position, magnetic moment) are set by the methods in `Class_unit_cell`.
-   `nn_conn`, `snd_nn_conn`: Enum values used to tag connection types.

## Dependencies/Interactions

-   **`Class_atom`**: Fundamental. The `unit_cell` type contains an array of `atom` objects. Methods in `Class_unit_cell` call methods of the `atom` type (e.g., `init_ferro_z`, `set_sphere`, `set_m_cart`, `free_atm`).
-   **`Class_helper`**: Provides utility functions like `my_norm2` (vector norm), `cross_prod`, `R_mtx` (rotation matrix generation), `n_times_phi` (for multi-q spirals/skyrmions), `deg_30`, `deg_60`, `PI`, `pos_eps`.
-   **`m_config`**: Used to read configuration parameters from a file (`CFG_t`, `CFG_get`, `CFG_get_size`).
-   **`output`**: For error messages (`error_msg`).
-   **`m_npy`**: Used for saving unit cell data to `.npy` files (`save_npy`).
-   **`mpi`**: Standard MPI module for parallel communication.
-   **`mypi`**: Custom MPI utilities, possibly for `MYPI_INT` or error checking (`check_ierr`).
-   **`Constants`**: Provides physical constants and mathematical constants (e.g., `PI`, `A_site`, `B_site` enumerators if site types are defined there).
-   **`class_Units`**: Provides the `units` type for managing physical units of quantities like length and energy.

This module acts as the primary store and manager of the crystal and magnetic structure, providing this information to other modules like `Class_hamiltonian` which then build the electronic model.
