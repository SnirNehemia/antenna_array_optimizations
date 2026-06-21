# ══════════════════════════════════════════════════════════════════
# COMPARE_CLASSICAL
# Compare optimised array weights against classical tapering and
# steering techniques for a synthetic N×N uniform rectangular array.
#
# Configuration is read from test_config.yaml at the project root.
# Two element-pattern sources are supported:
#   "synthetic" — ideal isotropic URA phase model (no real data needed)
#   "folder"    — CST Studio far-field exports loaded from a directory
#
# Results (figure + console log + config snapshot) are saved to a
# timestamped subfolder under results/compare_classical/.
#
# Usage:
#   python scripts/compare_classical.py                          # default config
#   python scripts/compare_classical.py --config path/cfg.yaml  # explicit config
#   python scripts/compare_classical.py --n 8                   # override n_side
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────
import argparse
import shutil
import sys
from datetime import datetime
from pathlib import Path

import warnings

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal.windows as spwin
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.cost.cost_function import angular_window_mask, compute_array_factor
from src.optimize.optimizer import run_optimizer

# ────────────────────────── BUILT-IN DEFAULTS ─────────────────────
# Applied for any key missing from the YAML config.

_DEFAULT_N_SIDE                = 4
_DEFAULT_D_OVER_LAMBDA         = 0.5
_DEFAULT_THETA_STEP_DEG        = 1.0
_DEFAULT_PHI_STEP_DEG          = 5.0
_DEFAULT_N_RESTARTS            = 5
_DEFAULT_MAX_ITERATIONS        = 500
_DEFAULT_COST_TOLERANCE        = 1e-9
_DEFAULT_PLOT_DYNAMIC_RANGE_DB = 60.0
_DEFAULT_ELEMENT_SOURCE        = "synthetic"
_DEFAULT_POLARIZATION          = "copol"

# Default scenarios — used when the config has no "scenarios" key.
# Scenario A: broadside peak, no steering.
# phi = 180°, phi_width = 360° ensures angular_window_mask covers all azimuth
# without wrapping artefacts at phi = 0°.
_DEFAULT_SCENARIO_A = {
    "label":           "Scenario A - Broadside SLL suppression",
    "steer_theta_deg": 0.0,
    "steer_phi_deg":   0.0,
    "null_theta_deg":  None,
    "directives": [
        {
            "type":        "peak",
            "theta":        7.5,
            "phi":          180.0,
            "theta_width":  15.0,
            "phi_width":    360.0,
            "weight":       1.0,
        },
    ],
}

# Scenario B: steered beam at θ = 20° with a null at θ = 60°.
_DEFAULT_SCENARIO_B = {
    "label":           "Scenario B - Steered 20deg, null at 60deg",
    "steer_theta_deg": 20.0,
    "steer_phi_deg":   0.0,
    "null_theta_deg":  60.0,
    "directives": [
        {
            "type":        "peak",
            "theta":        20.0,
            "phi":          0.0,
            "theta_width":  20.0,
            "phi_width":    30.0,
            "weight":       1.0,
        },
        {
            "type":        "null",
            "theta":        60.0,
            "phi":          0.0,
            "theta_width":  15.0,
            "phi_width":    20.0,
            "weight":       3.0,
        },
    ],
}

# Ordered list of classical technique labels (controls print/plot order).
CLASSICAL_TECHNIQUE_NAMES = [
    "uniform",
    "hamming",
    "hanning",
    "kaiser_b3",
    "kaiser_b6",
    "chebyshev_25dB",
    "chebyshev_40dB",
    "taylor_25dB",
]

# Colour / line style for each technique in the pattern overlay plot.
TECHNIQUE_PLOT_STYLE: dict = {
    "uniform":        {"color": "black",       "linestyle": "-",  "linewidth": 1.2, "alpha": 0.75},
    "hamming":        {"color": "royalblue",   "linestyle": "--", "linewidth": 1.2, "alpha": 0.85},
    "hanning":        {"color": "deepskyblue", "linestyle": ":",  "linewidth": 1.2, "alpha": 0.85},
    "kaiser_b3":      {"color": "seagreen",    "linestyle": "-.", "linewidth": 1.2, "alpha": 0.85},
    "kaiser_b6":      {"color": "limegreen",   "linestyle": "--", "linewidth": 1.2, "alpha": 0.85},
    "chebyshev_25dB": {"color": "darkorange",  "linestyle": "-",  "linewidth": 1.2, "alpha": 0.85},
    "chebyshev_40dB": {"color": "firebrick",   "linestyle": "--", "linewidth": 1.2, "alpha": 0.85},
    "taylor_25dB":    {"color": "purple",      "linestyle": ":",  "linewidth": 1.2, "alpha": 0.85},
    "optimized":      {"color": "magenta",     "linestyle": "-",  "linewidth": 2.5},
}


# ────────────────────────── CONFIG LOADING ────────────────────────

def _load_test_config(config_path):
    """Load test configuration from a YAML file, applying built-in defaults for missing keys.

    Args:
        config_path (Path): Path to the YAML configuration file.

    Returns:
        dict: Merged configuration with all required keys present.
    """
    cfg = {}
    if config_path.exists():
        with open(config_path, "r", encoding="utf-8") as fh:
            loaded = yaml.safe_load(fh) or {}
        cfg.update(loaded)
    else:
        print(
            f"WARNING: config file not found at '{config_path}'. "
            "Using built-in defaults."
        )

    cfg.setdefault("n_side",                _DEFAULT_N_SIDE)
    cfg.setdefault("d_over_lambda",         _DEFAULT_D_OVER_LAMBDA)
    cfg.setdefault("theta_step_deg",        _DEFAULT_THETA_STEP_DEG)
    cfg.setdefault("phi_step_deg",          _DEFAULT_PHI_STEP_DEG)
    cfg.setdefault("n_restarts",            _DEFAULT_N_RESTARTS)
    cfg.setdefault("max_iterations",        _DEFAULT_MAX_ITERATIONS)
    cfg.setdefault("cost_tolerance",        _DEFAULT_COST_TOLERANCE)
    cfg.setdefault("plot_dynamic_range_db", _DEFAULT_PLOT_DYNAMIC_RANGE_DB)
    cfg.setdefault("element_source",        _DEFAULT_ELEMENT_SOURCE)
    cfg.setdefault("element_patterns_dir",  "data/element_patterns/")
    cfg.setdefault("polarization",          _DEFAULT_POLARIZATION)
    cfg.setdefault("scenarios",             [_DEFAULT_SCENARIO_A, _DEFAULT_SCENARIO_B])

    return cfg


# ────────────────────────── ELEMENT PATTERN SYNTHESIS ─────────────

def build_ura_element_patterns(n_side, d_over_lambda, theta_deg, phi_deg):
    """Build synthetic element-pattern stack for an N×N uniform rectangular array.

    Each element is modelled as an isotropic radiator at position
    (m_c · d, n_c · d, 0) in the XY plane, where m_c and n_c are
    centre-referenced integer indices.  The element pattern for element
    (m, n) is the free-space progressive phase factor:

        E_{mn}(θ, φ) = exp(j · k · d · (m_c · sinθ · cosφ + n_c · sinθ · sinφ))

    With d = λ/2 this simplifies to exp(j · π · (m_c · u + n_c · v)),
    where u = sinθ · cosφ and v = sinθ · sinφ are the standard direction cosines.

    Args:
        n_side (int): Number of elements per side; total elements = n_side².
        d_over_lambda (float): Element spacing in wavelengths. Units: dimensionless.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.

    Returns:
        np.ndarray: Complex element patterns, shape (n_side², N_theta, N_phi).
            Units: V/m (1 V/m normalisation per element).
    """
    n_elements = n_side * n_side
    n_theta    = len(theta_deg)
    n_phi      = len(phi_deg)

    theta_rad = np.deg2rad(theta_deg)   # (N_theta,)
    phi_rad   = np.deg2rad(phi_deg)     # (N_phi,)

    # Direction cosines on the full (theta, phi) grid.
    sin_theta = np.sin(theta_rad)                                        # (N_theta,)
    cos_phi   = np.cos(phi_rad)                                          # (N_phi,)
    sin_phi   = np.sin(phi_rad)                                          # (N_phi,)
    u_grid    = sin_theta[:, np.newaxis] * cos_phi[np.newaxis, :]        # (N_theta, N_phi)
    v_grid    = sin_theta[:, np.newaxis] * sin_phi[np.newaxis, :]        # (N_theta, N_phi)
    # [MATLAB] [PHI, THETA] = meshgrid(phi_rad, theta_rad);
    # [MATLAB] u_grid = sin(THETA) .* cos(PHI);  v_grid = sin(THETA) .* sin(PHI);

    # Phase constant: k · d = 2π/λ · d = 2π · d_over_lambda.
    kd = 2.0 * np.pi * d_over_lambda

    # Centre-referenced element index offset so the array is symmetric about the origin.
    offset = (n_side - 1) / 2.0

    element_patterns_complex = np.empty((n_elements, n_theta, n_phi), dtype=complex)
    elem_idx = 0
    for m in range(n_side):
        for n in range(n_side):
            m_c = float(m) - offset     # centred x-index
            n_c = float(n) - offset     # centred y-index

            # Progressive phase for element at (m_c·d, n_c·d).
            phase_rad = kd * (m_c * u_grid + n_c * v_grid)    # (N_theta, N_phi)
            # [MATLAB] phase_rad = kd * (m_c * u_grid + n_c * v_grid);

            element_patterns_complex[elem_idx] = np.exp(1j * phase_rad)
            # [MATLAB] element_patterns_complex(elem_idx,:,:) = exp(1i * phase_rad);
            elem_idx += 1

    return element_patterns_complex


def _load_folder_element_patterns(cfg):
    """Load CST element patterns from a directory and stack them into a single array.

    Reads all ``.txt`` files from ``cfg["element_patterns_dir"]`` using the
    project CST parser, selects the field component specified by
    ``cfg["polarization"]``, and stacks the results into a (n_elements,
    N_theta, N_phi) array suitable for :func:`run_scenario`.

    The theta/phi grids are taken directly from the CST files.  The
    ``theta_step_deg`` / ``phi_step_deg`` config keys are ignored when
    element_source is "folder".

    Args:
        cfg (dict): Loaded test configuration (see :func:`_load_test_config`).

    Returns:
        tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray, int]:
            - ``element_patterns_complex``, shape (n_elements, N_theta, N_phi).
            - ``element_patterns_secondary``: cross-pol stack of the same shape when
              ``polarization == "total"``; ``None`` otherwise.
            - ``theta_deg``, shape (N_theta,). Units: degrees.
            - ``phi_deg``, shape (N_phi,). Units: degrees.
            - ``n_side`` (int): inferred from sqrt(n_elements).

    Raises:
        FileNotFoundError: If the directory contains no ``.txt`` files.
        ValueError: If the number of files is not a perfect square, or if
            an unsupported polarization is requested.
    """
    from src.io.cst_parser import load_element_patterns, get_component

    root         = Path(__file__).resolve().parents[1]
    patterns_dir = root / cfg["element_patterns_dir"]
    polarization = cfg.get("polarization", _DEFAULT_POLARIZATION)

    print(f"  Loading element patterns from: {patterns_dir}")
    raw_patterns = load_element_patterns(patterns_dir)
    n_elements   = len(raw_patterns)
    print(f"  Loaded {n_elements} element patterns.")

    # Infer array dimension.
    n_side_f = np.sqrt(n_elements)
    if not np.isclose(n_side_f, round(n_side_f)):
        raise ValueError(
            f"Number of CST element files ({n_elements}) is not a perfect square. "
            "compare_classical.py assumes an N x N square URA."
        )
    n_side = int(round(n_side_f))

    # Grid vectors come from the first file (the parser validates shape consistency).
    theta_deg = raw_patterns[0]["theta_deg"]    # (N_theta,)
    phi_deg   = raw_patterns[0]["phi_deg"]      # (N_phi,)

    # Select and stack the field component(s).
    # "total" loads both detected components so the cost function can use the
    # physically correct power sum |AF_a|² + |AF_b|² instead of the
    # approximate E_abs (magnitude-only) approach.
    component_names = sorted(raw_patterns[0]["components"].keys())
    ep_secondary = None
    if polarization.lower() == "total":
        if len(component_names) != 2:
            raise ValueError(
                f"polarization: 'total' requires exactly 2 polarization "
                f"components, found {component_names}."
            )
        name_primary, name_secondary = component_names
        ep          = np.stack([get_component(p, name_primary)["complex"]   for p in raw_patterns], axis=0)
        ep_secondary = np.stack([get_component(p, name_secondary)["complex"] for p in raw_patterns], axis=0)
        # Normalise both stacks together so the power ratio between components is preserved.
        peak_amp = float(np.max([np.max(np.abs(ep)), np.max(np.abs(ep_secondary))]))
        if peak_amp > 1e-30:
            ep          = ep          / peak_amp
            ep_secondary = ep_secondary / peak_amp
        return ep, ep_secondary, theta_deg, phi_deg, n_side
    else:
        matches = [name for name in component_names if name.lower() == polarization.lower()]
        if not matches:
            raise ValueError(
                f"Unknown polarization '{polarization}'. "
                f"Available components: {component_names} (or 'total')."
            )
        ep = np.stack([get_component(p, polarization)["complex"] for p in raw_patterns], axis=0)

    # Normalise to peak amplitude = 1 across all elements and angles.
    # CST exports have simulation-dependent absolute amplitudes; this makes the
    # dBi scale consistent: a coherent N-element array is bounded by 10·log10(N).
    peak_amp = float(np.max(np.abs(ep)))
    if peak_amp > 1e-30:
        ep = ep / peak_amp

    return ep, None, theta_deg, phi_deg, n_side


# ────────────────────────── STEERING PHASE ────────────────────────

def _steering_phase_vector(n_side, d_over_lambda, steer_theta_deg, steer_phi_deg):
    """Compute the 2D steering phase vector for a URA.

    Returns the complex conjugate of the array manifold vector evaluated at
    the desired steering direction (θ₀, φ₀).  Multiplying element weights by
    this vector shifts the main beam to (θ₀, φ₀).

    Args:
        n_side (int): Array dimension (N×N elements).
        d_over_lambda (float): Element spacing in wavelengths. Units: dimensionless.
        steer_theta_deg (float): Steering elevation angle. Units: degrees.
        steer_phi_deg (float): Steering azimuth angle. Units: degrees.

    Returns:
        np.ndarray: Complex steering phases, shape (n_side²,). Units: dimensionless.
    """
    steer_theta_rad = np.deg2rad(steer_theta_deg)
    steer_phi_rad   = np.deg2rad(steer_phi_deg)
    kd              = 2.0 * np.pi * d_over_lambda
    offset          = (n_side - 1) / 2.0

    # u₀ and v₀: direction cosines of the steering direction.
    u0 = np.sin(steer_theta_rad) * np.cos(steer_phi_rad)
    v0 = np.sin(steer_theta_rad) * np.sin(steer_phi_rad)
    # [MATLAB] u0 = sin(steer_theta_rad) * cos(steer_phi_rad);
    # [MATLAB] v0 = sin(steer_theta_rad) * sin(steer_phi_rad);

    steering_phases = np.empty(n_side * n_side, dtype=complex)
    elem_idx = 0
    for m in range(n_side):
        for n in range(n_side):
            m_c = float(m) - offset
            n_c = float(n) - offset

            # Conjugate of the manifold phase → cancels the phase delay at (θ₀, φ₀).
            phase_rad = kd * (m_c * u0 + n_c * v0)
            steering_phases[elem_idx] = np.exp(-1j * phase_rad)
            # [MATLAB] steering_phases(elem_idx) = exp(-1i * kd * (m_c * u0 + n_c * v0));
            elem_idx += 1

    return steering_phases


def _data_driven_steering_vector(element_patterns_complex, theta_deg, phi_deg,
                                  steer_theta_deg, steer_phi_deg):
    """Compute a model-free steering vector from actual element pattern data.

    Evaluates each element's far-field at the nearest grid point to the target
    direction and returns the conjugate unit-phase vector.  This is the correct
    approach for real antenna patterns: it works regardless of whether the CST
    exports encode geometric position phase relative to the array origin or to
    each element's own feed point.

    For synthetic (isotropic phase-only) element patterns this produces an
    identical result to :func:`_steering_phase_vector`.

    Args:
        element_patterns_complex (np.ndarray): Element patterns,
            shape (n_elements, N_theta, N_phi). Units: V/m.
        theta_deg (np.ndarray): Elevation grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth grid, shape (N_phi,). Units: degrees.
        steer_theta_deg (float): Steering elevation angle. Units: degrees.
        steer_phi_deg (float): Steering azimuth angle. Units: degrees.

    Returns:
        np.ndarray: Unit-amplitude complex steering phases, shape (n_elements,).
            Units: dimensionless.
    """
    theta_idx = int(np.argmin(np.abs(theta_deg - steer_theta_deg)))
    phi_idx   = int(np.argmin(np.abs(phi_deg   - steer_phi_deg)))

    # Array manifold: element patterns evaluated at the steering direction.
    manifold = element_patterns_complex[:, theta_idx, phi_idx]   # (n_elements,)
    # [MATLAB] manifold = element_patterns_complex(:, theta_idx, phi_idx);

    # Conjugate unit-phase: cancels each element's response phase at target direction.
    amp = np.maximum(np.abs(manifold), 1e-30)
    return np.conj(manifold) / amp
    # [MATLAB] steering = conj(manifold) ./ max(abs(manifold), 1e-30);


# ────────────────────────── CLASSICAL WEIGHTS ─────────────────────

def _window_1d(name, n_side):
    """Return a 1-D window of length n_side for the named technique.

    Args:
        name (str): Window name; one of the entries in CLASSICAL_TECHNIQUE_NAMES.
        n_side (int): Number of samples (= array elements per dimension).

    Returns:
        np.ndarray: Real window values, shape (n_side,). Units: dimensionless.

    Raises:
        ValueError: If name is not a recognised technique.
    """
    if name == "uniform":
        return np.ones(n_side, dtype=float)
    if name == "hamming":
        return np.hamming(n_side)
    if name == "hanning":
        return np.hanning(n_side)
    if name == "kaiser_b3":
        return np.kaiser(n_side, beta=3.0)
    if name == "kaiser_b6":
        return np.kaiser(n_side, beta=6.0)
    if name == "chebyshev_25dB":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            return spwin.chebwin(n_side, at=25.0)
    if name == "chebyshev_40dB":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            return spwin.chebwin(n_side, at=40.0)
    if name == "taylor_25dB":
        try:
            return spwin.taylor(n_side, nbar=4, sll=25, norm=True)
        except AttributeError:
            print(
                "WARNING: scipy.signal.windows.taylor requires scipy >= 1.7.0. "
                "Falling back to Chebyshev 25 dB for 'taylor_25dB'."
            )
            return spwin.chebwin(n_side, at=25.0)
    raise ValueError(
        f"Unknown classical technique '{name}'. "
        f"Valid names: {CLASSICAL_TECHNIQUE_NAMES}."
    )


def classical_weights(name, n_side, d_over_lambda,
                      steer_theta_deg=0.0, steer_phi_deg=0.0,
                      element_patterns_complex=None, theta_deg=None, phi_deg=None):
    """Compute complex element weights for a classical tapering + steering method.

    The 2-D aperture taper is the separable outer product of two identical 1-D
    windows applied along the x and y dimensions of the array.  A progressive
    steering phase is then applied element-wise to point the main beam toward
    (steer_theta_deg, steer_phi_deg).

    Steering mode is selected automatically:

    - **Data-driven** (element_patterns_complex, theta_deg, phi_deg all provided):
      evaluates the actual array manifold at the target direction and conjugates
      its phase.  Works correctly for real CST element patterns regardless of
      whether the exports encode geometric position phase relative to the array
      origin or to each element's own feed point.  Produces identical results to
      the geometric model for synthetic (phase-only) patterns.
    - **Geometric model** (element patterns not provided): computes progressive
      phase shifts from the assumed element positions and d_over_lambda.  Only
      accurate for synthetic patterns or when the CST exports are phased relative
      to the array centre.

    Args:
        name (str): Classical technique name (see CLASSICAL_TECHNIQUE_NAMES).
        n_side (int): Array dimension; total elements = n_side².
        d_over_lambda (float): Element spacing in wavelengths. Units: dimensionless.
        steer_theta_deg (float): Steering elevation angle. Units: degrees.
        steer_phi_deg (float): Steering azimuth angle. Units: degrees.
        element_patterns_complex (np.ndarray | None): Element patterns,
            shape (n_elements, N_theta, N_phi). Triggers data-driven steering
            when not None. Units: V/m.
        theta_deg (np.ndarray | None): Elevation grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray | None): Azimuth grid, shape (N_phi,). Units: degrees.

    Returns:
        np.ndarray: Complex element weights, shape (n_side²,). Units: dimensionless.
    """
    window_1d = _window_1d(name, n_side)                                          # (n_side,)
    window_2d = window_1d[:, np.newaxis] * window_1d[np.newaxis, :]               # (n_side, n_side)
    # [MATLAB] window_2d = window_1d * window_1d';

    weights_taper_complex = window_2d.ravel().astype(complex)                     # (n_side²,)
    # [MATLAB] weights_taper_complex = complex(window_2d(:), zeros(n_side^2, 1));

    use_data_driven = (
        element_patterns_complex is not None
        and theta_deg is not None
        and phi_deg is not None
    )
    if use_data_driven:
        steering_phases = _data_driven_steering_vector(
            element_patterns_complex, theta_deg, phi_deg, steer_theta_deg, steer_phi_deg
        )
    else:
        steering_phases = _steering_phase_vector(
            n_side, d_over_lambda, steer_theta_deg, steer_phi_deg
        )

    return weights_taper_complex * steering_phases                                 # (n_side²,)
    # [MATLAB] weights_complex = weights_taper_complex .* steering_phases;


# ────────────────────────── PATTERN CUTS ──────────────────────────

def principal_plane_cut(array_factor_complex, phi_deg):
    """Extract the principal plane cut from a 2-D complex array factor.

    Combines the φ = 0° half-plane (positive-θ side) and the φ = 180°
    half-plane (negative-θ side) into a single cut spanning −90° to +90°.
    The display convention is:
        theta_display_deg < 0  →  comes from the φ = 180° half-plane
        theta_display_deg > 0  →  comes from the φ = 0° half-plane

    Args:
        array_factor_complex (np.ndarray): Complex array factor,
            shape (N_theta, N_phi). Units: V/m.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.

    Returns:
        np.ndarray: Linear magnitude of the principal-plane cut,
            shape (2·N_theta − 1,). Units: V/m.
    """
    phi0_idx   = int(np.argmin(np.abs(phi_deg - 0.0)))
    phi180_idx = int(np.argmin(np.abs(phi_deg - 180.0)))

    af_pos_vm = np.abs(array_factor_complex[:, phi0_idx])    # (N_theta,)  phi=0°
    af_neg_vm = np.abs(array_factor_complex[:, phi180_idx])  # (N_theta,)  phi=180°

    cut_vm = np.concatenate([af_neg_vm[::-1][:-1], af_pos_vm])   # (2·N_theta − 1,)
    # [MATLAB] cut_vm = [af_neg_vm(end:-1:2); af_pos_vm];
    return cut_vm


def principal_plane_theta_axis(theta_deg):
    """Build the symmetric theta display axis for a principal-plane cut.

    Mirrors the one-sided theta grid [0°, 90°] to produce a two-sided
    axis [−90°, ..., 0°, ..., +90°] of length 2·N_theta − 1.

    Args:
        theta_deg (np.ndarray): One-sided elevation grid, shape (N_theta,).
            Must start at 0°. Units: degrees.

    Returns:
        np.ndarray: Full display axis, shape (2·N_theta − 1,). Units: degrees.
    """
    # [MATLAB] theta_display_deg = [-fliplr(theta_deg(2:end)), theta_deg];
    return np.concatenate([-theta_deg[::-1][:-1], theta_deg])   # (2·N_theta − 1,)


# ────────────────────────── METRICS ───────────────────────────────

def _power_normalize_weights(weights_complex):
    """Normalise complex weights so their total power equals 1.

    Divides by the L2 norm of the weight vector, matching the internal
    normalisation applied by build_cost_function.  After normalisation the
    peak array-factor magnitude equals sqrt(N) for a fully coherent uniform
    array, giving 10·log10(N) dBi — the theoretical maximum directivity for
    N isotropic elements.

    Args:
        weights_complex (np.ndarray): Raw complex weights, shape (N,).
            Units: dimensionless.

    Returns:
        np.ndarray: Power-normalised weights, shape (N,). Units: dimensionless.
    """
    total_power = float(np.sum(np.abs(weights_complex) ** 2, axis=0))
    safe_power  = max(total_power, 1e-30)
    return weights_complex / np.sqrt(safe_power)
    # [MATLAB] weights_complex = weights_complex / sqrt(max(sum(abs(weights_complex).^2), 1e-30));


def evaluate_directive_metrics(af_norm_mag_vm, theta_deg, phi_deg, directive):
    """Evaluate min / solid-angle-weighted-mean / max pattern level in one directive window.

    The angular window is taken from the directive dict — the same region the
    optimizer uses — so all techniques are evaluated on identical ground.
    The ``peak_dbi`` field reflects the true beamforming gain because the
    input AF must already be computed from power-normalised weights
    (see ``_power_normalize_weights``).

    Args:
        af_norm_mag_vm (np.ndarray): |AF| from power-normalised weights,
            shape (N_theta, N_phi). Units: V/m (1 V/m per element after norm).
        theta_deg (np.ndarray): Elevation grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth grid, shape (N_phi,). Units: degrees.
        directive (dict): Directive dict with keys ``type``, ``theta``,
            and either ``width`` or ``theta_width`` + ``phi_width``.
            Optional ``phi`` (default 0.0). Units: degrees.

    Returns:
        dict: Keys:
            - ``directive_type`` (str): ``"peak"`` or ``"null"``.
            - ``peak_dbi`` (float): Global |AF| peak in dBi after power
              normalisation.  For uniform weights = 10·log10(N). Units: dBi.
            - ``min_db`` (float): Minimum |AF| in window, relative to peak.
              Units: dB.
            - ``mean_db`` (float): Solid-angle-weighted mean |AF| in window,
              relative to peak.  Units: dB.
            - ``max_db`` (float): Maximum |AF| in window, relative to peak.
              Units: dB.
    """
    directive_type = directive["type"]
    target_theta   = directive["theta"]
    target_phi     = directive.get("phi", 0.0)
    sym_width      = directive.get("width", None)
    theta_width    = directive.get("theta_width", sym_width)
    phi_width      = directive.get("phi_width",   sym_width)

    mask = angular_window_mask(
        theta_deg, phi_deg, target_theta, target_phi, theta_width, phi_width,
    )
    # [MATLAB] mask = angular_window_mask(theta_deg, phi_deg, target_theta, target_phi, theta_width, phi_width);

    global_peak_vm = float(np.max(af_norm_mag_vm))
    peak_dbi       = 20.0 * np.log10(max(global_peak_vm, 1e-30))

    window_values_vm = af_norm_mag_vm[mask]   # (N_in_window,)

    if window_values_vm.size == 0:
        nan = float("nan")
        return {"directive_type": directive_type, "peak_dbi": peak_dbi,
                "min_db": nan, "mean_db": nan, "max_db": nan}

    min_vm = float(np.min(window_values_vm))
    max_vm = float(np.max(window_values_vm))
    min_db = 20.0 * np.log10(max(min_vm, 1e-30)) - peak_dbi
    max_db = 20.0 * np.log10(max(max_vm, 1e-30)) - peak_dbi
    # [MATLAB] min_db = 20*log10(max(min(window_values_vm), 1e-30)) - peak_dbi;

    # Solid-angle-weighted mean |AF| in the window.
    sin_theta       = np.sin(np.deg2rad(theta_deg))                # (N_theta,)
    weights_sa      = mask * sin_theta[:, np.newaxis]              # (N_theta, N_phi)
    total_weight_sa = float(np.sum(weights_sa))
    mean_vm         = float(np.sum(af_norm_mag_vm * weights_sa)) / max(total_weight_sa, 1e-30)
    mean_db         = 20.0 * np.log10(max(mean_vm, 1e-30)) - peak_dbi
    # [MATLAB] mean_vm = sum(af_norm_mag_vm(:) .* weights_sa(:)) / max(sum(weights_sa(:)), 1e-30);

    return {
        "directive_type": directive_type,
        "peak_dbi":       peak_dbi,
        "min_db":         min_db,
        "mean_db":        mean_db,
        "max_db":         max_db,
    }


# ────────────────────────── SCENARIO RUNNER ───────────────────────

def run_scenario(scenario, element_patterns_complex, theta_deg, phi_deg,
                 n_side, d_over_lambda, cfg,
                 element_patterns_secondary=None):
    """Compute classical and optimised patterns for one scenario.

    Args:
        scenario (dict): Scenario definition dict.
        element_patterns_complex (np.ndarray): URA element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m.
        theta_deg (np.ndarray): Elevation grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth grid, shape (N_phi,). Units: degrees.
        n_side (int): Array dimension.
        d_over_lambda (float): Element spacing in wavelengths. Units: dimensionless.
        cfg (dict): Full test configuration (supplies optimizer settings).
        element_patterns_secondary (np.ndarray | None): Cross-pol element patterns,
            shape (N_elements, N_theta, N_phi). When provided, pattern power is
            computed as ``|AF_copol|² + |AF_cross|²``. Default: None.

    Returns:
        dict: Mapping technique_name → {
            ``weights_complex``   (np.ndarray, shape (N_elements,)),
            ``cut_vm``            (np.ndarray, shape (N_cut,)),
            ``directive_metrics`` (list[dict], one entry per directive),
        }.
    """
    steer_theta_deg = scenario["steer_theta_deg"]
    steer_phi_deg   = scenario["steer_phi_deg"]
    directives      = scenario["directives"]

    results: dict = {}

    def _process(weights_complex):
        """Compute AF magnitude (or total-power magnitude), principal-plane cut, and per-directive metrics."""
        weights_norm = _power_normalize_weights(weights_complex)
        af_copol     = compute_array_factor(weights_norm, element_patterns_complex)
        if element_patterns_secondary is not None:
            af_cross  = compute_array_factor(weights_norm, element_patterns_secondary)
            # Effective magnitude: sqrt(|AF_copol|² + |AF_cross|²).
            # np.abs() on this real non-negative array is a no-op, so
            # principal_plane_cut and evaluate_directive_metrics work unchanged.
            af_mag_vm = np.sqrt(np.abs(af_copol) ** 2 + np.abs(af_cross) ** 2)
            # [MATLAB] af_mag_vm = sqrt(abs(af_copol).^2 + abs(af_cross).^2);
        else:
            af_mag_vm = np.abs(af_copol)
        cut_vm      = principal_plane_cut(af_mag_vm, phi_deg)   # (N_cut,)
        dir_metrics = [
            evaluate_directive_metrics(af_mag_vm, theta_deg, phi_deg, d)
            for d in directives
        ]
        return weights_complex, cut_vm, dir_metrics

    # ── Classical techniques ───────────────────────────────────────
    for name in CLASSICAL_TECHNIQUE_NAMES:
        weights_raw = classical_weights(
            name, n_side, d_over_lambda, steer_theta_deg, steer_phi_deg,
            element_patterns_complex, theta_deg, phi_deg,
        )
        weights_complex, cut_vm, dir_metrics = _process(weights_raw)
        results[name] = {
            "weights_complex":   weights_complex,
            "cut_vm":            cut_vm,
            "directive_metrics": dir_metrics,
        }

    # ── Optimised technique ────────────────────────────────────────
    n_restarts = cfg.get("n_restarts", _DEFAULT_N_RESTARTS)
    optimizer_config = {
        "max_iterations":          cfg.get("max_iterations",  _DEFAULT_MAX_ITERATIONS),
        "cost_tolerance":          cfg.get("cost_tolerance",  _DEFAULT_COST_TOLERANCE),
        "n_restarts":              n_restarts,
        "use_uniform_init":        True,
        "use_single_element_init": False,
    }
    print(f"  Running optimizer ({n_restarts} restarts)...")
    opt_result = run_optimizer(
        element_patterns_complex, theta_deg, phi_deg, directives, optimizer_config,
        element_patterns_secondary=element_patterns_secondary,
    )
    weights_opt, cut_opt_vm, dir_metrics_opt = _process(opt_result["weights_complex"])
    results["optimized"] = {
        "weights_complex":   weights_opt,
        "cut_vm":            cut_opt_vm,
        "directive_metrics": dir_metrics_opt,
    }

    return results


# ────────────────────────── CONSOLE OUTPUT ────────────────────────

def print_metrics_table(scenario_label, results, directives):
    """Print per-directive min / mean / max tables for all techniques.

    For peak directives the Min column is marked (*) — it is the weakest
    point inside the intended coverage area.  For null directives the Max
    column is marked (*) — it is the worst-case leakage inside the null region.

    Args:
        scenario_label (str): Human-readable scenario name.
        results (dict): Output of run_scenario.
        directives (list[dict]): Scenario directives list (from the scenario dict).
    """
    technique_order = list(CLASSICAL_TECHNIQUE_NAMES) + ["optimized"]
    techniques      = [n for n in technique_order if n in results]

    print(f"\n{'=' * 75}")
    print(f"  {scenario_label}")
    print(f"{'=' * 75}")

    for d_idx, directive in enumerate(directives):
        d_type   = directive["type"].upper()
        t_center = directive["theta"]
        p_center = directive.get("phi", 0.0)
        sym_w    = directive.get("width", None)
        t_width  = directive.get("theta_width", sym_w)

        print(f"\n  Directive {d_idx + 1}: {d_type}  "
              f"theta={t_center:.0f}deg  phi={p_center:.0f}deg  "
              f"theta_width={t_width:.0f}deg")

        if directive["type"] == "peak":
            h_min, h_max = "Min* (dB)", "Max (dB)"
            note = "  * Min = worst-case coverage point in window"
        else:
            h_min, h_max = "Min (dB)", "Max* (dB)"
            note = "  * Max = worst-case leakage point in window"

        print(f"  {'Technique':<18} {'Peak (dBi)':>11} {h_min:>10} {'Mean (dB)':>10} {h_max:>10}")
        print(f"  {'-'*18} {'-'*11} {'-'*10} {'-'*10} {'-'*10}")

        for name in techniques:
            dm     = results[name]["directive_metrics"][d_idx]
            peak_s = f"{dm['peak_dbi']:>11.1f}"
            mean_s = f"{dm['mean_db']:>10.1f}" if not np.isnan(dm["mean_db"]) else f"{'N/A':>10}"

            if not np.isnan(dm["min_db"]):
                min_s = f"{dm['min_db']:>9.1f}*" if directive["type"] == "peak" else f"{dm['min_db']:>10.1f}"
            else:
                min_s = f"{'N/A':>10}"

            if not np.isnan(dm["max_db"]):
                max_s = f"{dm['max_db']:>9.1f}*" if directive["type"] == "null" else f"{dm['max_db']:>10.1f}"
            else:
                max_s = f"{'N/A':>10}"

            marker = " <--" if name == "optimized" else ""
            print(f"  {name:<18}{peak_s}{min_s}{mean_s}{max_s}{marker}")

        print(note)


# ────────────────────────── PLOTTING ──────────────────────────────

def _plot_pattern_overlay(ax, scenario, results, theta_display_deg, dynamic_range_db):
    """Draw the absolute-gain principal-plane cut overlay for one scenario.

    The y-axis shows actual gain in dBi (power-normalised array factor, so
    uniform weights peak at 10·log10(N) dBi for isotropic elements).  All
    traces share the same dBi scale so gain loss from tapering is immediately
    visible.

    Args:
        ax (matplotlib.axes.Axes): Target axes.
        scenario (dict): Scenario definition dict.
        results (dict): Output of run_scenario.
        theta_display_deg (np.ndarray): Symmetric theta display axis. Units: degrees.
        dynamic_range_db (float): Floor below the global peak. Units: dB.
    """
    null_theta_deg  = scenario["null_theta_deg"]
    technique_order = list(CLASSICAL_TECHNIQUE_NAMES) + ["optimized"]

    # First pass: convert to dBi and find the global peak for y-axis range.
    cuts_dbi     = {}
    max_peak_dbi = -np.inf
    for name in technique_order:
        if name not in results:
            continue
        cut_vm = results[name]["cut_vm"]
        with np.errstate(divide="ignore", invalid="ignore"):
            cut_dbi = 20.0 * np.log10(np.maximum(cut_vm, 1e-30))
        cuts_dbi[name] = cut_dbi
        max_peak_dbi   = max(max_peak_dbi, float(np.max(cut_dbi)))

    floor_dbi = max_peak_dbi - dynamic_range_db

    # Second pass: plot clipped traces.
    for name in technique_order:
        if name not in cuts_dbi:
            continue
        cut_plot = np.maximum(cuts_dbi[name], floor_dbi)
        style    = TECHNIQUE_PLOT_STYLE.get(name, {})
        ax.plot(theta_display_deg, cut_plot, label=name, **style)

    # -3 dB reference line relative to the global peak (uniform = max for isotropic).
    ax.axhline(max_peak_dbi - 3.0, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)

    if null_theta_deg is not None:
        ax.axvline(
            null_theta_deg,
            color="red", linestyle=":", linewidth=1.2,
            label=f"null target {null_theta_deg:.0f} deg",
        )

    ax.set_xlim([float(theta_display_deg[0]), float(theta_display_deg[-1])])
    ax.set_ylim([floor_dbi, max_peak_dbi + 2.0])
    ax.set_xlabel("theta display (degrees)")
    ax.set_ylabel("Gain (dBi)")
    ax.set_title(scenario["label"] + " - pattern cut (phi = 0 deg plane)")
    ax.legend(fontsize=7, loc="lower center", ncol=2, framealpha=0.85)
    ax.grid(True, alpha=0.25)


def _plot_directive_whiskers(ax, directive, directive_idx, results):
    """Draw a min / mean / max whisker chart for one directive across all techniques.

    Each technique is one point on the x-axis.  The dot shows the
    solid-angle-weighted mean inside the directive window; the vertical whisker
    spans from min to max.  The critical endpoint is highlighted:

    - Peak directive: the **min** tip is marked (worst-case coverage point).
    - Null  directive: the **max** tip is marked (worst-case leakage point).

    The optimised technique is drawn in magenta; classical ones in steelblue.

    Args:
        ax (matplotlib.axes.Axes): Target axes.
        directive (dict): One directive dict from the scenario.
        directive_idx (int): Index into each technique's directive_metrics list.
        results (dict): Output of run_scenario.
    """
    technique_order = list(CLASSICAL_TECHNIQUE_NAMES) + ["optimized"]
    names           = [n for n in technique_order if n in results]
    d_type          = directive["type"]

    x         = np.arange(len(names), dtype=float)
    peak_dbis = [results[n]["directive_metrics"][directive_idx]["peak_dbi"] for n in names]
    means_dbi = [results[n]["directive_metrics"][directive_idx]["mean_db"] + peak_dbis[i]
                 for i, n in enumerate(names)]
    mins_dbi  = [results[n]["directive_metrics"][directive_idx]["min_db"]  + peak_dbis[i]
                 for i, n in enumerate(names)]
    maxs_dbi  = [results[n]["directive_metrics"][directive_idx]["max_db"]  + peak_dbis[i]
                 for i, n in enumerate(names)]

    for i, name in enumerate(names):
        is_opt = (name == "optimized")
        stem_c = "magenta"  if is_opt else "steelblue"
        lw     = 2.2        if is_opt else 1.4

        ax.plot([x[i], x[i]], [mins_dbi[i], maxs_dbi[i]],
                color=stem_c, linewidth=lw, alpha=0.85, zorder=3)
        ax.scatter([x[i]], [means_dbi[i]], marker="o", s=45,
                   color=stem_c, zorder=5, linewidths=0)

    # Highlighted critical tip: min for peak directive, max for null directive.
    crit_vals, crit_marker, crit_label = (
        (mins_dbi, "v", "Min (worst coverage) *") if d_type == "peak"
        else (maxs_dbi, "^", "Max (worst leakage) *")
    )

    text_va     = "top"    if d_type == "peak" else "bottom"
    text_offset = -5       if d_type == "peak" else 5        # points
    for i, name in enumerate(names):
        tip_c = "crimson" if (name == "optimized") else "darkorange"
        ax.scatter([x[i]], [crit_vals[i]], marker=crit_marker, s=80,
                   color=tip_c, zorder=6, linewidths=0)
        ax.annotate(
            f"{crit_vals[i]:.1f}",
            xy=(x[i], crit_vals[i]),
            xytext=(0, text_offset),
            textcoords="offset points",
            ha="center", va=text_va,
            fontsize=6, color=tip_c, zorder=7,
        )

    t_center = directive["theta"]
    p_center = directive.get("phi", 0.0)
    sym_w    = directive.get("width", None)
    t_width  = directive.get("theta_width", sym_w)
    ax.set_title(
        f"{d_type.upper()} window: "
        f"theta={t_center:.0f}deg  phi={p_center:.0f}deg  width={t_width:.0f}deg",
        fontsize=8,
    )
    ax.set_ylabel("Gain (dBi)", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=35, ha="right", fontsize=7)
    ax.axhline(max(peak_dbis), color="black", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.grid(True, axis="y", alpha=0.25)

    from matplotlib.lines import Line2D
    ax.legend(
        handles=[
            Line2D([0], [0], color="steelblue", linewidth=1.4, label="[min, max] range"),
            Line2D([0], [0], marker="o", color="steelblue", linestyle="None",
                   markersize=6, label="Mean (SA-weighted)"),
            Line2D([0], [0], marker=crit_marker, color="darkorange", linestyle="None",
                   markersize=8, label=crit_label),
            Line2D([0], [0], color="magenta", linewidth=2.2, label="optimized"),
        ],
        fontsize=7, loc="best",
    )


def plot_comparison(all_scenario_results, theta_deg, n_side, cfg, output_path):
    """Render the full comparison figure (N_scenarios × 2 columns each).

    Left column: absolute-gain (dBi) principal-plane pattern overlay.
    Right column: one whisker chart per directive, stacked vertically —
        dot = solid-angle-weighted mean, whisker = [min, max],
        critical tip highlighted (min for peak, max for null).

    Args:
        all_scenario_results (list[tuple[dict, dict]]): List of (scenario, results)
            pairs, one per scenario.
        theta_deg (np.ndarray): One-sided elevation grid. Units: degrees.
        n_side (int): Array dimension, used in figure title.
        cfg (dict): Full test configuration (supplies d_over_lambda, element_source,
            plot_dynamic_range_db).
        output_path (Path): File path to save the figure.
    """
    theta_display_deg = principal_plane_theta_axis(theta_deg)
    d_over_lambda     = cfg["d_over_lambda"]
    dynamic_range_db  = cfg["plot_dynamic_range_db"]
    elem_label        = "CST elements" if cfg["element_source"] == "folder" else "isotropic elements"

    max_dir = max(len(s["directives"]) for s, _ in all_scenario_results)
    row_h   = max(4.5 * max_dir, 5.5)
    n_rows  = len(all_scenario_results)
    fig     = plt.figure(figsize=(16, row_h * n_rows))
    fig.suptitle(
        f"Classical tapering / steering vs. L-BFGS-B optimiser  -  "
        f"{n_side}x{n_side} URA  (d = {d_over_lambda} lambda, {elem_label})",
        fontsize=12, fontweight="bold",
    )

    outer_gs = gridspec.GridSpec(n_rows, 2, figure=fig, hspace=0.5, wspace=0.35)

    for row_idx, (scenario, results) in enumerate(all_scenario_results):
        ax_pat = fig.add_subplot(outer_gs[row_idx, 0])
        _plot_pattern_overlay(ax_pat, scenario, results, theta_display_deg, dynamic_range_db)

        n_dir    = len(scenario["directives"])
        inner_gs = gridspec.GridSpecFromSubplotSpec(
            n_dir, 1, subplot_spec=outer_gs[row_idx, 1], hspace=0.55,
        )
        for d_idx, directive in enumerate(scenario["directives"]):
            ax_dir = fig.add_subplot(inner_gs[d_idx])
            _plot_directive_whiskers(ax_dir, directive, d_idx, results)

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved -> {output_path}")
    plt.show()


# ────────────────────────── MAIN ──────────────────────────────────

class _Tee:
    """Forward write() calls to both the original stdout and a log file."""

    def __init__(self, real_stdout, log_file):
        self._stdout  = real_stdout
        self._logfile = log_file

    def write(self, text):
        self._stdout.write(text)
        self._logfile.write(text)

    def flush(self):
        self._stdout.flush()
        self._logfile.flush()


def _create_output_dir(base_dir):
    """Create and return a timestamped subfolder under base_dir."""
    timestamp  = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    output_dir = Path(base_dir) / timestamp
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _parse_args():
    """Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed argument namespace.
    """
    default_cfg = Path(__file__).resolve().parents[1] / "test_config.yaml"
    parser = argparse.ArgumentParser(
        description=(
            "Compare classical tapering / steering techniques against "
            "the L-BFGS-B optimiser on an N x N URA. "
            "All settings are read from test_config.yaml."
        )
    )
    parser.add_argument(
        "--config", type=Path, default=default_cfg,
        metavar="PATH",
        help=f"Path to the test configuration YAML file. Default: {default_cfg}.",
    )
    parser.add_argument(
        "--n", type=int, default=None,
        metavar="N",
        help="Override n_side from the config (N x N elements).",
    )
    return parser.parse_args()


def main():
    """Entry point: load config, build the array, run scenarios, print and plot results."""
    args = _parse_args()
    cfg  = _load_test_config(args.config)

    # --n overrides config (backward compatibility).
    if args.n is not None:
        cfg["n_side"] = args.n

    # ── Create timestamped output directory ────────────────────────
    project_root = Path(__file__).resolve().parents[1]
    output_dir   = _create_output_dir(project_root / "results" / "compare_classical")
    print(f"Output directory: {output_dir}")

    # ── Tee stdout → log file for the rest of the run ─────────────
    log_path        = output_dir / "console_output.txt"
    original_stdout = sys.stdout

    with open(log_path, "w", encoding="utf-8") as log_file:
        sys.stdout = _Tee(original_stdout, log_file)
        try:
            d_over_lambda  = cfg["d_over_lambda"]
            element_source = cfg["element_source"]
            polarization   = cfg.get("polarization", _DEFAULT_POLARIZATION)

            # ── Element patterns ───────────────────────────────────
            if element_source == "folder":
                element_patterns_complex, element_patterns_secondary, theta_deg, phi_deg, n_side = \
                    _load_folder_element_patterns(cfg)
                cfg["n_side"] = n_side    # keep cfg consistent for figure title
                n_elements    = element_patterns_complex.shape[0]
                print(
                    f"\nLoaded {n_side}x{n_side} URA from folder "
                    f"({n_elements} elements, {len(theta_deg)} theta x {len(phi_deg)} phi points)."
                )
            elif element_source == "synthetic":
                if polarization.lower() == "total":
                    raise ValueError(
                        "polarization: 'total' is not supported with element_source: 'synthetic'. "
                        "Synthetic patterns have no cross-polarization component. "
                        "Use element_source: 'folder' with real CST data, or set "
                        "polarization: 'copol'."
                    )
                n_side    = cfg["n_side"]
                theta_deg = np.arange(0.0, 90.0 + cfg["theta_step_deg"], cfg["theta_step_deg"])
                phi_deg   = np.arange(0.0, 360.0, cfg["phi_step_deg"])
                n_elements = n_side * n_side
                print(
                    f"\nBuilding {n_side}x{n_side} URA element patterns "
                    f"({n_elements} elements, d = {d_over_lambda} lambda, "
                    f"{len(theta_deg)} theta x {len(phi_deg)} phi points) ..."
                )
                element_patterns_complex = build_ura_element_patterns(
                    n_side, d_over_lambda, theta_deg, phi_deg
                )
                element_patterns_secondary = None
                print(f"  Element patterns shape: {element_patterns_complex.shape}")
            else:
                raise ValueError(
                    f"Unknown element_source '{element_source}'. "
                    "Valid options: 'synthetic', 'folder'."
                )

            # ── Scenario loop ──────────────────────────────────────
            n_side   = cfg["n_side"]
            all_scenario_results = []
            for scenario in cfg["scenarios"]:
                print(f"\n{'=' * 68}")
                print(f"  {scenario['label']}")
                print(f"{'=' * 68}")
                results = run_scenario(
                    scenario, element_patterns_complex,
                    theta_deg, phi_deg, n_side, d_over_lambda, cfg,
                    element_patterns_secondary=element_patterns_secondary,
                )
                print_metrics_table(scenario["label"], results, scenario["directives"])
                all_scenario_results.append((scenario, results))

            # ── Figure ────────────────────────────────────────────
            figure_path = output_dir / f"compare_classical_{n_side}x{n_side}.png"
            plot_comparison(all_scenario_results, theta_deg, n_side, cfg, figure_path)

            # ── Config snapshot ────────────────────────────────────
            shutil.copy(args.config, output_dir / "test_config.yaml")
            print(f"\nFiles saved to: {output_dir}")

        finally:
            sys.stdout = original_stdout


if __name__ == "__main__":
    main()
