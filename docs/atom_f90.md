# Documentation for `src/atom.f90`

## Overview

The `Class_atom` module defines a derived data type `atom` which represents a single atom within a crystal unit cell. It stores information about the atom's position, its magnetic moment (orientation), its site type (e.g., for bipartite lattices like honeycomb), and details about its neighboring atoms, including their indices and the vectors connecting to them. This module is fundamental for constructing the unit cell and subsequently the Hamiltonian of the system.

## Key Components

### `module Class_atom`

This is the main module.

### `enum, bind(c)` for Site Types

Defines enumerators for classifying atom sites, particularly useful for lattices with sublattices:
-   `A_site`: Represents an atom on the A sublattice.
-   `B_site`: Represents an atom on the B sublattice.
-   `no_site`: Used when sublattice distinction is not applicable or needed.

### `type atom`

A derived data type to store properties of an individual atom.

**Key Members (Parameters):**

-   `m_phi (real(8))`: Azimuthal angle $\phi$ of the atom's magnetic moment in spherical coordinates.
-   `m_theta (real(8))`: Polar angle $\theta$ of the atom's magnetic moment in spherical coordinates. The convention usually follows physics standards (theta is angle from z-axis, phi is angle from x-axis in xy-plane).
-   `pos (real(8), dimension(3))`: Cartesian coordinates $(x, y, z)$ of the atom's position in real space, typically in atomic units or units defined by a lattice constant.
-   `site_type (integer(4))`: An integer representing the site type, using values from the `enum` (A_site, B_site, no_site).
-   `neigh_idx (integer, allocatable :: neigh_idx(:))`: An array storing the indices (within the unit cell's atom list) of its neighboring atoms.
-   `neigh_conn (real(8), allocatable :: neigh_conn(:,:))`: A 2D array where each row is a 3D vector representing the real-space connection vector from this atom to a corresponding neighbor in `neigh_idx`.
-   `conn_type (integer(4), allocatable :: conn_type(:))`: An array storing the type of connection for each neighbor (e.g., nearest neighbor, second nearest neighbor, from an `enum` likely defined in `Class_unit_cell`).
-   `me (integer)`: MPI rank of the current process.
-   `nProcs (integer)`: Total number of MPI processes.

**Key Procedures (Methods):**

-   **`init_ferro_z(p_pos, site) result(self)`**: Constructor function. Initializes an `atom` object.
    -   **Arguments:**
        -   `p_pos (real(8), intent(in))`: The 3D position of the atom.
        -   `site (integer, optional)`: The site type of the atom. Defaults to `no_site` if not provided.
    -   Sets the initial magnetic moment to point along the positive z-axis (phi=0, theta=0).
    -   Assigns the given position and site type.
    -   Initializes MPI rank and process count.
-   **`get_m_cart(self) result(coord)`**: Function that returns the Cartesian coordinates $(m_x, m_y, m_z)$ of the atom's magnetic moment, calculated from its spherical angles `m_phi` and `m_theta`. Assumes a unit magnetic moment.
-   **`set_sphere(self, phi, theta)`**: Subroutine to set the atom's magnetic moment using spherical coordinates.
    -   **Arguments:**
        -   `phi (real(8), intent(in))`: The new azimuthal angle.
        -   `theta (real(8), intent(in))`: The new polar angle.
-   **`set_m_cart(self, x, y, z)`**: Subroutine to set the atom's magnetic moment using Cartesian coordinates.
    -   **Arguments:**
        -   `x, y, z (real(8), intent(in))`: The Cartesian components of the new magnetic moment.
    -   Converts these Cartesian coordinates to spherical angles (`m_phi`, `m_theta`) and stores them.
    -   Includes a check to ensure the input vector is approximately normalized.
-   **`free_atm(self)`**: Subroutine to deallocate allocatable members of the `atom` object (`neigh_idx`, `neigh_conn`). Note: `conn_type` is missing from deallocation in the provided code, which could be a potential memory leak if it's allocated.
-   **`compare_to_root(self) result(success)`**: Function used for testing in an MPI environment. It compares all members of the `atom` object on the current process with the corresponding values from the root process (rank 0) and returns `.True.` if they all match within tolerance, `.False.` otherwise. This helps ensure data consistency across MPI processes.

## Important Variables/Constants (within `type atom`)

-   **`m_phi, m_theta (real(8))`**: These two angles define the orientation of the atom's magnetic moment, which is crucial for magnetic Hamiltonians.
-   **`pos(3) (real(8))`**: The atom's position, essential for calculating distances and phase factors in hopping terms.
-   **`site_type (integer(4))`**: Distinguishes atoms on different sublattices, which can have different onsite energies or behave differently under certain interactions.
-   **`neigh_idx(:)`**: Array of indices pointing to neighboring atoms. This, along with `neigh_conn` and `conn_type`, defines the local environment and connectivity of the atom, which is necessary for constructing hopping terms in the Hamiltonian.

## Dependencies/Interactions

-   **`Class_helper`**:
    -   `my_norm2`: Used in `set_m_cart` to check normalization of the input magnetic vector and in `compare_to_root` to compare position vectors.
    -   `mtx_norm`: Used in `compare_to_root` to compare neighbor connection matrices.
    -   `error_msg`: Used for reporting errors, for example, if the spin vector in `set_m_cart` is not normalized, or if inconsistencies are found in `compare_to_root`.
-   **`Constants`**: While not explicitly shown using a constant from this module in the provided snippet, it's listed as a dependency. It likely provides physical constants or mathematical constants (like PI if not intrinsic) that might be used in more extensive calculations involving atoms. The site type enum itself could also be considered part of a constants collection.
-   **`mpi`**: Standard MPI module used for:
    -   Getting MPI communicator size (`MPI_Comm_size`) and rank (`MPI_Comm_rank`) in `init_ferro_z`.
    -   Broadcasting data (`MPI_Bcast`) and performing comparisons in `compare_to_root`.
    -   Error checking for MPI calls (`check_ierr`, which is likely a wrapper from `Class_helper` or `mypi`).

The `atom` type is a building block for the `unit_cell` type in `Class_unit_cell.f90`. The `unit_cell` module will typically create an array of `atom` objects and then populate their neighbor information based on the lattice geometry.
