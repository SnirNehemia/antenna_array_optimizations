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

# Fixed leading columns shared by every CST far-field export.
COL_THETA = 0  # Elevation angle [deg]
COL_PHI   = 1  # Azimuth angle [deg]

# Regex matching one "field" column header: Abs(<name>), Phase(<name>), or Ax.Ratio.
# Matches are scanned left-to-right starting at column index 2 (after Theta, Phi).
_FIELD_HEADER_PATTERN = re.compile(r"Abs\(([^)]*)\)|Phase\(([^)]*)\)|Ax\.Ratio")


# ────────────────────────── HELPER FUNCTIONS ──────────────────────

def _parse_header_columns(header_line):
    """Identify the role of each data column from the CST header line.

    Searches the header for ``Abs(<name>)`` / ``Phase(<name>)`` pairs (e.g.
    ``Abs(Copol)``/``Phase(Copol)``, or ``Abs(Theta)``/``Phase(Theta)``) and the
    ``Ax.Ratio`` column. Columns 0 and 1 are always ``Theta`` and ``Phi``. The
    first ``Abs(<name>)`` column with no matching ``Phase(<name>)`` is treated
    as the magnitude-only total-field column (``Abs(E)`` or ``Abs(Grlz)``).

    Args:
        header_line (str): First line of the CST export file (column labels).

    Returns:
        dict: With keys:

            - ``e_abs_col`` (int): Column index of the magnitude-only total
              field column.
            - ``axial_ratio_col`` (int): Column index of ``Ax.Ratio``.
            - ``components`` (dict[str, dict[str, int]]): Maps each detected
              polarization-component name (as written in the header, e.g.
              ``"Copol"``, ``"Cross"``, ``"Theta"``, ``"Phi"``) to
              ``{"abs_col": int, "phase_col": int}``.

    Raises:
        ValueError: If the header does not contain exactly the expected set
            of columns (one magnitude-only column, one ``Ax.Ratio`` column,
            and at least one complete Abs/Phase pair).
    """
    abs_cols = {}    # name -> column index
    phase_cols = {}  # name -> column index
    e_abs_col = None
    axial_ratio_col = None

    col_idx = COL_PHI + 1  # field columns start right after Theta, Phi
    for match in _FIELD_HEADER_PATTERN.finditer(header_line):
        if match.group(0) == "Ax.Ratio":
            axial_ratio_col = col_idx
        elif match.group(1) is not None:
            name = match.group(1).strip()
            if name in abs_cols:
                raise ValueError(
                    f"Duplicate 'Abs({name})' column in header: {header_line!r}"
                )
            abs_cols[name] = col_idx
        else:
            name = match.group(2).strip()
            if name in phase_cols:
                raise ValueError(
                    f"Duplicate 'Phase({name})' column in header: {header_line!r}"
                )
            phase_cols[name] = col_idx
        col_idx += 1

    # The magnitude-only total-field column is the Abs(...) entry with no
    # matching Phase(...) entry (e.g. Abs(E), Abs(Grlz)).
    e_abs_names = [name for name in abs_cols if name not in phase_cols]
    if len(e_abs_names) != 1:
        raise ValueError(
            f"Expected exactly one magnitude-only Abs(...) column (e.g. "
            f"'Abs(E)' or 'Abs(Grlz)'), found {e_abs_names} in header: "
            f"{header_line!r}"
        )
    e_abs_col = abs_cols.pop(e_abs_names[0])

    if axial_ratio_col is None:
        raise ValueError(f"No 'Ax.Ratio' column found in header: {header_line!r}")

    components = {}
    for name, abs_col in abs_cols.items():
        if name not in phase_cols:
            raise ValueError(
                f"'Abs({name})' column has no matching 'Phase({name})' column "
                f"in header: {header_line!r}"
            )
        components[name] = {"abs_col": abs_col, "phase_col": phase_cols[name]}

    if not components:
        raise ValueError(
            f"No Abs(<name>)/Phase(<name>) component pairs found in header: "
            f"{header_line!r}"
        )

    return {
        "e_abs_col": e_abs_col,
        "axial_ratio_col": axial_ratio_col,
        "components": components,
    }


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


def get_component(pattern, polarization_name):
    """Look up a polarization component's complex grid by name (case-insensitive).

    Args:
        pattern (dict): Element pattern dict returned by :func:`parse_cst_file`.
        polarization_name (str): Component name to look up, e.g. ``"Copol"``,
            ``"Cross"``, ``"Theta"``, ``"Phi"``. Matched case-insensitively
            against ``pattern["components"]``.

    Returns:
        dict: ``{"abs": ..., "phase": ..., "complex": ...}`` for the matched
        component (same arrays as stored in ``pattern["components"]``).

    Raises:
        KeyError: If no component matches ``polarization_name``, listing the
            available component names.
    """
    target = polarization_name.strip().lower()
    for name, comp in pattern["components"].items():
        if name.lower() == target:
            return comp
    available = sorted(pattern["components"].keys())
    raise KeyError(
        f"Polarization component '{polarization_name}' not found. "
        f"Available components in this data: {available}."
    )


# ────────────────────────── MAIN FUNCTIONS ────────────────────────

def parse_cst_file(filepath):
    """Parse a single CST Studio far-field ASCII export file.

    Reads an 8-column tab/space-separated `.txt` file produced by CST Studio,
    auto-detects the angular grid resolution, and returns all field components
    reshaped into (N_theta, N_phi) grids. Polarization components are detected
    generically from the header's ``Abs(<name>)``/``Phase(<name>)`` column
    pairs (e.g. ``Copol``/``Cross``, or ``Theta``/``Phi``), so this function
    supports any CST far-field export basis without code changes.

    For each detected component ``<name>``, the complex field is built as:

        complex(theta, phi) = Abs(<name>) * exp(j * Phase(<name>) [rad])

    Args:
        filepath (str | Path): Path to the CST far-field export `.txt` file.

    Returns:
        dict: Parsed element pattern with keys:

            - ``theta_deg``    (np.ndarray, shape (N_theta,)) — elevation angle grid.
              Units: degrees.
            - ``phi_deg``      (np.ndarray, shape (N_phi,)) — azimuth angle grid.
              Units: degrees.
            - ``E_abs``        (np.ndarray, shape (N_theta, N_phi)) — magnitude-only
              total-field column (``Abs(E)`` or ``Abs(Grlz)``). Units: as exported.
            - ``axial_ratio``  (np.ndarray, shape (N_theta, N_phi)) — axial ratio.
              Units: dimensionless.
            - ``components``   (dict[str, dict]) — one entry per detected
              ``Abs(<name>)``/``Phase(<name>)`` pair, keyed by ``<name>`` exactly
              as written in the header (e.g. ``"Copol"``, ``"Cross"``,
              ``"Theta"``, ``"Phi"``). Each value is a dict with:

                - ``abs``     (np.ndarray, shape (N_theta, N_phi)) — magnitude.
                - ``phase``   (np.ndarray, shape (N_theta, N_phi)) — phase, degrees.
                - ``complex`` (np.ndarray, shape (N_theta, N_phi), complex) —
                  ``abs * exp(j * phase_rad)``.

    Raises:
        FileNotFoundError: If ``filepath`` does not exist.
        ValueError: If the file has an unexpected number of columns, the
            header columns cannot be identified, or the angular grid is
            incomplete.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"CST export file not found: {filepath}")

    with open(filepath, "r", encoding="utf-8") as fh:
        header_line = fh.readline()

    header_info = _parse_header_columns(header_line)

    # Load all numeric data, skipping the 2-line header.
    raw_data = np.loadtxt(filepath, skiprows=N_HEADER_LINES)
    # [MATLAB] raw_data = dlmread(filepath, '', N_HEADER_LINES, 0);

    if raw_data.ndim != 2 or raw_data.shape[1] != N_DATA_COLUMNS:
        raise ValueError(
            f"Expected {N_DATA_COLUMNS} columns in '{filepath.name}', "
            f"got shape {raw_data.shape}."
        )

    # Extract flat angle columns.
    theta_deg_flat = raw_data[:, COL_THETA]
    phi_deg_flat   = raw_data[:, COL_PHI]

    n_theta, n_phi = _detect_grid_shape(theta_deg_flat, phi_deg_flat)

    # Extract sorted 1-D angle vectors from the unique values.
    theta_deg = np.unique(theta_deg_flat)   # shape (N_theta,), sorted ascending
    phi_deg   = np.unique(phi_deg_flat)     # shape (N_phi,),   sorted ascending
    # [MATLAB] theta_deg = unique(theta_deg_flat);
    # [MATLAB] phi_deg   = unique(phi_deg_flat);

    def _reshape(col_idx):
        # Theta is the fast (inner) axis: each block of N_theta consecutive rows
        # shares the same Phi value. Reshape to (N_phi, N_theta) first, then
        # transpose to arrive at the canonical (N_theta, N_phi) layout.
        # [MATLAB] reshape(flat, n_theta, n_phi)
        return raw_data[:, col_idx].reshape(n_phi, n_theta).T

    e_abs_grid       = _reshape(header_info["e_abs_col"])
    axial_ratio_grid = _reshape(header_info["axial_ratio_col"])

    # Build the complex field for every detected polarization component.
    # Formula: complex = Abs(<name>) * exp(j * Phase(<name>) [rad])
    # [MATLAB] phase_rad_grid = phase_deg_grid .* (pi / 180);
    # [MATLAB] complex_grid   = abs_grid .* exp(1j .* phase_rad_grid);
    components = {}
    for name, cols in header_info["components"].items():
        abs_grid       = _reshape(cols["abs_col"])
        phase_deg_grid = _reshape(cols["phase_col"])
        phase_rad_grid = phase_deg_grid * (np.pi / 180.0)
        complex_grid   = abs_grid * np.exp(1j * phase_rad_grid)
        components[name] = {
            "abs":     abs_grid,
            "phase":   phase_deg_grid,
            "complex": complex_grid,
        }

    return {
        "theta_deg":   theta_deg,
        "phi_deg":     phi_deg,
        "E_abs":       e_abs_grid,
        "axial_ratio": axial_ratio_grid,
        "components":  components,
    }


def load_element_patterns(directory_path):
    """Load and parse all CST element pattern files from a directory.

    Discovers every ``.txt`` file in the given directory, sorts them by the
    element index encoded in the filename (e.g. ``[1]``, ``[2]``, …), and
    parses each one with :func:`parse_cst_file`. All elements must share the
    same angular grid shape and the same set of polarization components.

    Args:
        directory_path (str | Path): Path to the directory containing CST
            export ``.txt`` files.

    Returns:
        list[dict]: One dict per element (see :func:`parse_cst_file` for the
            dict schema), sorted by element index in ascending order.

    Raises:
        FileNotFoundError: If ``directory_path`` contains no ``.txt`` files.
        ValueError: If any file cannot have its element index extracted, if
            the elements do not all share the same (N_theta, N_phi) grid
            shape, or if they do not all expose the same polarization
            components.
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
    reference_components = None

    for txt_file in txt_files:
        pattern = parse_cst_file(txt_file)

        grid_shape = pattern["E_abs"].shape
        component_names = set(pattern["components"].keys())
        if reference_shape is None:
            reference_shape = grid_shape
            reference_components = component_names
        elif grid_shape != reference_shape:
            raise ValueError(
                f"Grid shape mismatch across elements: '{txt_file.name}' has shape "
                f"{grid_shape}, expected {reference_shape}. All elements must share "
                "the same (N_theta, N_phi) angular grid."
            )
        elif component_names != reference_components:
            raise ValueError(
                f"Polarization component mismatch across elements: '{txt_file.name}' "
                f"has components {sorted(component_names)}, expected "
                f"{sorted(reference_components)}. All elements must share the same "
                "set of polarization components."
            )

        element_patterns.append(pattern)

    return element_patterns
