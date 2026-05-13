# ══════════════════════════════════════════════════════════════════
# CST_PARSER
# Parse CST Studio far-field export .txt files into structured numpy arrays.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────
import re
from pathlib import Path

import numpy as np

# ────────────────────────── CONSTANTS ─────────────────────────────

# Number of header lines to skip at the top of every CST export file.
# Line 1: column labels. Line 2: dash separator.
N_HEADER_LINES = 2

# Number of data columns per row in the CST export.
N_DATA_COLUMNS = 8

# Column indices matching the CST far-field ASCII export format.
COL_THETA       = 0  # Elevation angle [deg]
COL_PHI         = 1  # Azimuth angle [deg]
COL_E_ABS       = 2  # Total E-field magnitude [V/m]
COL_CROSS_ABS   = 3  # Cross-polarization magnitude [V/m]
COL_CROSS_PHASE = 4  # Cross-polarization phase [deg]
COL_COPOL_ABS   = 5  # Co-polarization magnitude [V/m]
COL_COPOL_PHASE = 6  # Co-polarization phase [deg]
COL_AXIAL_RATIO = 7  # Axial ratio [dimensionless]


# ────────────────────────── HELPER FUNCTIONS ──────────────────────

def _detect_grid_shape(theta_deg_flat, phi_deg_flat):
    """Infer the (N_theta, N_phi) grid dimensions from flat angle arrays.

    Counts unique Theta and Phi values and verifies their product equals
    the total number of samples, confirming a complete rectangular grid.

    Args:
        theta_deg_flat (np.ndarray): Flat array of Theta values, shape (N_rows,).
            Units: degrees.
        phi_deg_flat (np.ndarray): Flat array of Phi values, shape (N_rows,).
            Units: degrees.

    Returns:
        tuple[int, int]: (n_theta, n_phi) — number of unique Theta and Phi values.

    Raises:
        ValueError: If n_theta * n_phi does not equal the total row count,
            indicating an incomplete or malformed angular grid.
    """
    n_theta = len(np.unique(theta_deg_flat))
    # [MATLAB] n_theta = length(unique(theta_deg_flat));
    n_phi = len(np.unique(phi_deg_flat))
    # [MATLAB] n_phi = length(unique(phi_deg_flat));

    n_rows = len(theta_deg_flat)
    if n_theta * n_phi != n_rows:
        raise ValueError(
            f"Grid shape mismatch: {n_theta} unique Theta × {n_phi} unique Phi "
            f"= {n_theta * n_phi}, but file has {n_rows} data rows. "
            "Expected a complete rectangular angular grid."
        )

    return n_theta, n_phi


def _extract_element_index(filepath):
    """Extract the element index integer from a CST filename.

    Matches the bracket notation used by CST Studio exports, e.g.
    ``farfield (f=1.4e+4) [3]_E_field.txt`` → 3.

    Args:
        filepath (str | Path): Path to the CST export file.

    Returns:
        int: Element index extracted from the filename.

    Raises:
        ValueError: If no ``[N]`` bracket pattern is found in the filename.
    """
    filename = Path(filepath).name
    match = re.search(r'\[(\d+)\]', filename)
    # [MATLAB] regexp(filename, '\[(\d+)\]', 'tokens') — returns cell of matches
    if match is None:
        raise ValueError(
            f"Could not extract element index from filename '{filename}'. "
            "Expected a pattern like '[N]' (e.g. '[3]') in the filename."
        )
    return int(match.group(1))


# ────────────────────────── MAIN FUNCTIONS ────────────────────────

def parse_cst_file(filepath):
    """Parse a single CST Studio far-field ASCII export file.

    Reads an 8-column tab/space-separated `.txt` file produced by CST Studio,
    auto-detects the angular grid resolution, and returns all field components
    reshaped into (N_theta, N_phi) grids. The complex element pattern is built
    from the co-polarization magnitude and phase as:

        E_complex(theta, phi) = Abs(Copol) * exp(j * Phase(Copol) [rad])

    Args:
        filepath (str | Path): Path to the CST far-field export `.txt` file.

    Returns:
        dict: Parsed element pattern with keys:

            - ``theta_deg``   (np.ndarray, shape (N_theta,)) — elevation angle grid.
              Units: degrees.
            - ``phi_deg``     (np.ndarray, shape (N_phi,)) — azimuth angle grid.
              Units: degrees.
            - ``E_abs``       (np.ndarray, shape (N_theta, N_phi)) — total E-field
              magnitude. Units: V/m.
            - ``copol_abs``   (np.ndarray, shape (N_theta, N_phi)) — co-pol magnitude.
              Units: V/m.
            - ``copol_phase`` (np.ndarray, shape (N_theta, N_phi)) — co-pol phase.
              Units: degrees.
            - ``E_complex``   (np.ndarray, shape (N_theta, N_phi), complex) — complex
              co-pol element pattern. Units: V/m.

    Raises:
        FileNotFoundError: If ``filepath`` does not exist.
        ValueError: If the file has an unexpected number of columns or an
            incomplete angular grid.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"CST export file not found: {filepath}")

    # Load all numeric data, skipping the 2-line header.
    raw_data = np.loadtxt(filepath, skiprows=N_HEADER_LINES)
    # [MATLAB] raw_data = dlmread(filepath, '', N_HEADER_LINES, 0);

    if raw_data.ndim != 2 or raw_data.shape[1] != N_DATA_COLUMNS:
        raise ValueError(
            f"Expected {N_DATA_COLUMNS} columns in '{filepath.name}', "
            f"got shape {raw_data.shape}."
        )

    # Extract flat angle and field columns.
    theta_deg_flat      = raw_data[:, COL_THETA]
    phi_deg_flat        = raw_data[:, COL_PHI]
    e_abs_flat          = raw_data[:, COL_E_ABS]
    copol_abs_flat      = raw_data[:, COL_COPOL_ABS]
    copol_phase_deg_flat = raw_data[:, COL_COPOL_PHASE]

    n_theta, n_phi = _detect_grid_shape(theta_deg_flat, phi_deg_flat)

    # Extract sorted 1-D angle vectors from the unique values.
    theta_deg = np.unique(theta_deg_flat)   # shape (N_theta,), sorted ascending
    phi_deg   = np.unique(phi_deg_flat)     # shape (N_phi,),   sorted ascending
    # [MATLAB] theta_deg = unique(theta_deg_flat);
    # [MATLAB] phi_deg   = unique(phi_deg_flat);

    # Reshape flat columns into (N_theta, N_phi) grids.
    # Theta is the fast (inner) axis: each block of N_theta consecutive rows
    # shares the same Phi value. Reshape to (N_phi, N_theta) first, then
    # transpose to arrive at the canonical (N_theta, N_phi) layout.
    e_abs_grid           = e_abs_flat.reshape(n_phi, n_theta).T
    copol_abs_grid       = copol_abs_flat.reshape(n_phi, n_theta).T
    copol_phase_deg_grid = copol_phase_deg_flat.reshape(n_phi, n_theta).T
    # [MATLAB] e_abs_grid           = reshape(e_abs_flat,           n_theta, n_phi);
    # [MATLAB] copol_abs_grid       = reshape(copol_abs_flat,       n_theta, n_phi);
    # [MATLAB] copol_phase_deg_grid = reshape(copol_phase_deg_flat, n_theta, n_phi);

    # Build the complex element pattern from co-polarization magnitude and phase.
    # Convert phase from degrees to radians for all internal complex arithmetic.
    # Formula: E_complex = |E_copol| * exp(j * phi_copol)
    copol_phase_rad_grid = copol_phase_deg_grid * (np.pi / 180.0)
    e_complex_grid = copol_abs_grid * np.exp(1j * copol_phase_rad_grid)
    # [MATLAB] copol_phase_rad_grid = copol_phase_deg_grid .* (pi / 180);
    # [MATLAB] e_complex_grid = copol_abs_grid .* exp(1j .* copol_phase_rad_grid);

    return {
        "theta_deg":   theta_deg,
        "phi_deg":     phi_deg,
        "E_abs":       e_abs_grid,
        "copol_abs":   copol_abs_grid,
        "copol_phase": copol_phase_deg_grid,
        "E_complex":   e_complex_grid,
    }


def load_element_patterns(directory_path):
    """Load and parse all CST element pattern files from a directory.

    Discovers every ``.txt`` file in the given directory, sorts them by the
    element index encoded in the filename (e.g. ``[1]``, ``[2]``, …), and
    parses each one with :func:`parse_cst_file`. All elements must share the
    same angular grid shape.

    Args:
        directory_path (str | Path): Path to the directory containing CST
            export ``.txt`` files.

    Returns:
        list[dict]: One dict per element (see :func:`parse_cst_file` for the
            dict schema), sorted by element index in ascending order.

    Raises:
        FileNotFoundError: If ``directory_path`` contains no ``.txt`` files.
        ValueError: If any file cannot have its element index extracted, or if
            the elements do not all share the same (N_theta, N_phi) grid shape.
    """
    directory_path = Path(directory_path)
    txt_files = sorted(
        directory_path.glob("*.txt"),
        key=_extract_element_index,
    )
    # [MATLAB] avoid list comprehension; use dir() + a for-loop to collect file paths

    if not txt_files:
        raise FileNotFoundError(
            f"No .txt files found in directory: {directory_path}"
        )

    element_patterns = []
    reference_shape = None

    for txt_file in txt_files:
        pattern = parse_cst_file(txt_file)

        grid_shape = pattern["E_complex"].shape
        if reference_shape is None:
            reference_shape = grid_shape
        elif grid_shape != reference_shape:
            raise ValueError(
                f"Grid shape mismatch across elements: '{txt_file.name}' has shape "
                f"{grid_shape}, expected {reference_shape}. All elements must share "
                "the same (N_theta, N_phi) angular grid."
            )

        element_patterns.append(pattern)

    return element_patterns
