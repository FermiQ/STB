# Documentation for `src/main.f90`

## Overview

This program, named `STB`, is the main driver for a solid-state physics calculation package. It initializes MPI for parallel processing, reads input configuration files, and then processes each file to perform various calculations like band structure, density of states (DOS), Berry curvature related quantities (e.g., Hall conductivity, orbital magnetization), and Anomalous Hall Conductivity (ACA). The program is designed to handle multiple input files and distribute calculations across MPI ranks.

## Key Components

### `program STB`

The main program unit.
- Initializes and finalizes MPI.
- Seeds the random number generator.
- Calls `get_inp_files` to determine the input configuration files to be processed.
- Broadcasts the number of files to all MPI ranks.
- Loops through each input file and calls `process_file` to handle the calculations for that specific configuration.

### `subroutine process_file(inp_file)`

This subroutine orchestrates the calculations for a single input configuration file.
- **Arguments:**
    - `inp_file (character(len=300), intent(in))`: The path to the input configuration file.
- Reads the configuration file using `CFG_read_file` (only on the root MPI rank).
- Calls `add_full_cfg` to populate the configuration with default or derived values.
- Retrieves various calculation flags (e.g., `perform_band`, `perform_dos`, `calc_hall`) from the configuration.
- Broadcasts these flags and other essential parameters (like `fermi_type`, `pert_log`) to all MPI ranks.
- Initializes the k-space representation (`Ksp`) using `init_k_space` based on the configuration.
- Saves the complete configuration to a file (e.g., `setup.cfg`) using `save_cfg` (only on the root MPI rank).
- Performs calculations based on the enabled flags:
    - Band structure (`Ksp%calc_and_print_band()`).
    - Sets Fermi energy if `fermi_type` is "fixed".
    - Density of States (`Ksp%calc_and_print_dos()`).
    - Finds Fermi energy based on filling if `fermi_type` is "filling" and DOS was calculated.
    - Berry quantities (Hall conductivity, orbital magnetization) (`Ksp%calc_berry_quantities(pert_log)`).
    - Anomalous Hall Conductivity (`Ksp%calc_ACA()`).
    - Plots Berry curvature (`Ksp%plot_omega_square()`) if `plot_omega` is true.
- Measures and prints execution time for initialization and total processing.
- Frees allocated memory associated with `Ksp` using `Ksp%free_ksp()`.
- Clears the configuration object `cfg` (only on the root MPI rank).

### `subroutine get_inp_files(n_files, inp_files)`

Determines the list of input configuration files to be processed based on command-line arguments.
- **Output Arguments:**
    - `n_files (integer, intent(out))`: The total number of input files.
    - `inp_files (character(len=300), allocatable, intent(out))`: An array of strings, where each string is the path to an input file.
- Logic for handling command-line arguments (only on the root MPI rank):
    - If 0 arguments: Prints an error and aborts.
    - If 1 argument: Assumes it's a single input file.
    - If 2 arguments: Assumes the first is a base string and the second is the number of files (`N`). Generates filenames like `base_str0.cfg`, `base_str1.cfg`, ..., `base_str(N-1).cfg`.
    - If 3 arguments: Assumes the first is a base string, the second is a start index, and the third is an end index. Generates filenames accordingly.
- Broadcasts `n_files` to all MPI ranks.
- Allocates `inp_files` on non-root ranks.
- Broadcasts each input filename to all MPI ranks.

### `subroutine add_full_cfg(cfg)`

Populates the configuration object (`cfg`) with a comprehensive set of default parameters and their descriptions for various aspects of the calculation (units, Hamiltonian parameters, grid details, band structure, DOS, Berry calculations, output settings, ACA, layer dropout, and plotting).
- **Arguments:**
    - `cfg (type(CFG_t))`: The configuration object to be populated.
- Uses `CFG_add` to define numerous configuration parameters with their default values and descriptions. This acts as a schema for the expected configuration options.

### `subroutine save_cfg(cfg)`

Saves the current configuration to a file.
- **Arguments:**
    - `cfg (type(CFG_t))`: The configuration object to be saved.
- Adds a timestamp for when the calculation started (`calculation%start_time`) to the configuration.
- Retrieves the output prefix (`output%band_prefix`) from the configuration.
- Writes the configuration to a file named `[prefix]setup.cfg` using `CFG_write`.

## Important Variables/Constants

- **`me (integer)`**: The MPI rank of the current process.
- **`root (integer, parameter, implicit in mpi module)`**: Typically rank 0, used for tasks like I/O and broadcasting initial data. The actual value is likely defined in the `mpi` module (e.g., `parameter :: root = 0`).
- **`Ksp (type(k_space))`**: An instance of the `k_space` class, which likely encapsulates the k-space grid, Hamiltonian, and methods for performing various calculations.
- **`cfg (type(CFG_t))`**: A derived type variable holding the configuration parameters for the current calculation.
- **`inp_files (character(len=300), allocatable :: inp_files(:))`**: Array storing the names of input configuration files.
- **`n_files (integer)`**: The number of input files to process.
- **`time_fmt (character(len=*), parameter)`**: Format string for timing output.

## Dependencies/Interactions

The `STB` program and its subroutines utilize several modules:

- **`Class_k_space`**: Provides the `k_space` derived type (`Ksp`) and associated methods for calculations (band structure, DOS, Berry quantities, ACA).
- **`m_config`**: Provides the `CFG_t` type and subroutines for managing configuration files (`CFG_read_file`, `CFG_get`, `CFG_add`, `CFG_write`, `CFG_clear`).
- **`m_npy`**: Likely used for reading/writing data in the NumPy .npy format, possibly for compatibility with Python-based tools or for efficient data storage. (Though not explicitly used in `main.f90` directly, it might be a dependency of `Class_k_space` or other modules).
- **`output`**: Potentially contains subroutines for formatted output or logging (e.g., `error_msg`).
- **`mpi`**: Provides MPI constants (e.g., `MPI_COMM_WORLD`, `MYPI_INT`, `MPI_LOGICAL`, `MPI_CHARACTER`) and procedures (`MPI_Init`, `MPI_Comm_rank`, `MPI_Bcast`, `MPI_Finalize`, `MPI_Wtime`). The constant `root` is also likely defined here.
- **`Constants`**: Might define physical constants or other shared parameters.
- **`Class_unit_cell`**: Likely provides the definition for the unit cell structure, which would be a component of the Hamiltonian within `Class_k_space`. (e.g. `Ksp%ham%UC%num_atoms`)
- **`mypi`**: Possibly a custom MPI utility module, or it could be defining `MYPI_INT` if it's not a standard MPI constant.

The program interacts with the file system by reading input configuration files and writing output files (e.g., band structure data, DOS data, and the `setup.cfg` file). It also interacts with the command line to receive input file names or patterns.
