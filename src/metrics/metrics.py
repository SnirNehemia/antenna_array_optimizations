# ══════════════════════════════════════════════════════════════════
# METRICS
# Post-optimization evaluation: per-directive gain scores and cost breakdown.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────
import numpy as np

from src.cost.cost_function import compute_array_factor, angular_window_mask

# ────────────────────────── CONSTANTS ─────────────────────────────

# Small offset added inside log10 to prevent log(0) for zero-power grid points.
LOG10_EPSILON = 1e-30

# Default azimuth target when not specified by a directive.
DEFAULT_TARGET_PHI_DEG = 0.0

# Default directive weight when not specified by the user.
DEFAULT_DIRECTIVE_WEIGHT = 1.0


# ────────────────────────── HELPERS ───────────────────────────────

def _nearest_index(grid_deg, target_deg):
    """Return the index of the nearest grid point to a target angle.

    Args:
        grid_deg (np.ndarray): 1D angle grid in degrees.
        target_deg (float): Target angle in degrees.

    Returns:
        int: Index of the grid point closest to target_deg.
    """
    return int(np.argmin(np.abs(grid_deg - target_deg)))
    # [MATLAB] [~, idx] = min(abs(grid_deg - target_deg));



# ────────────────────────── HELPERS (continued) ───────────────────

def _compute_directivity_dbi_grid(array_factor, theta_deg, phi_deg):
    """Compute directivity in dBi over the full (N_theta × N_phi) angular grid.

    Uses the discrete rectangle-rule approximation to the spherical integral:

        P_total ≈ Σ_θ Σ_φ |AF(θ,φ)|² · sin(θ) · Δθ · Δφ   [(V/m)² · sr]
        D(θ,φ)  =  4π · |AF(θ,φ)|² / P_total               [dimensionless]

    The absolute V/m² scale of the element patterns cancels exactly in the
    ratio, so the result is independent of CST normalisation.

    Args:
        array_factor (np.ndarray): Complex array factor, shape (N_theta, N_phi).
            Units: V/m (arbitrary absolute scale).
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
            Must cover 0° – 180° for a full-sphere integral.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.

    Returns:
        np.ndarray: Directivity grid, shape (N_theta, N_phi). Units: dBi.
    """
    power_grid = np.abs(array_factor) ** 2
    # [MATLAB] power_grid = abs(array_factor) .^ 2;

    theta_rad  = np.deg2rad(theta_deg)
    dtheta_rad = float(np.diff(theta_rad).mean())  if len(theta_rad) > 1 else np.pi
    dphi_rad   = float(np.deg2rad(np.diff(phi_deg).mean())) if len(phi_deg) > 1 else 2.0 * np.pi
    # [MATLAB] dtheta_rad = mean(diff(theta_rad)); dphi_rad = mean(diff(phi_rad));

    sin_theta = np.sin(theta_rad)  # shape (N_theta,); solid-angle weighting

    # Numerical spherical integral: Σ |AF|² sin(θ) Δθ Δφ
    p_total = np.sum(power_grid * sin_theta[:, np.newaxis]) * dtheta_rad * dphi_rad
    # [MATLAB] p_total = sum(sum(power_grid .* sin(theta_rad)')) * dtheta_rad * dphi_rad;

    directivity_linear = (4.0 * np.pi * power_grid) / p_total
    # [MATLAB] directivity_linear = 4 * pi * power_grid / p_total;

    return 10.0 * np.log10(directivity_linear + LOG10_EPSILON)
    # [MATLAB] directivity_dbi_grid = 10 * log10(directivity_linear + eps);


# ────────────────────────── PUBLIC INTERFACE ──────────────────────

def evaluate_metrics(element_patterns_stacked, theta_deg, phi_deg,
                     weights_complex, directives, cost_history):
    """Evaluate post-optimization performance metrics for all beam-shaping directives.

    Computes the array factor from the final optimized weights, then for each
    directive reports the realized power level at the target angle (in dB) and
    the cost contribution. Null depth is expressed relative to the global pattern peak.

    Args:
        element_patterns_stacked (np.ndarray): Complex element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        weights_complex (np.ndarray): Optimized complex weights, shape (N_elements,).
            Units: dimensionless (V/V).
        directives (list[dict]): Beam-shaping directives. Each dict must have:
            - ``type`` (str): ``"peak"`` or ``"null"``.
            - ``theta`` (float): Target elevation angle. Units: degrees.
            - ``width`` (float): Full angular window width. Units: degrees.
            Optional keys:
            - ``phi`` (float): Target azimuth angle (default 0.0). Units: degrees.
            - ``weight`` (float): Directive weight lambda_k (default 1.0).
        cost_history (list[float]): Cost value per iteration from the optimizer.

    Returns:
        dict: Metrics summary with keys:

            - ``total_cost`` (float): Composite cost J(x*) at the optimized weights.
            - ``iteration_count`` (int): Number of optimizer iterations recorded.
            - ``global_peak_dbi`` (float): Peak directivity over the full grid. Units: dBi.
            - ``directive_metrics`` (list[dict]): Per-directive results. Each entry has:

                - ``type`` (str): ``"peak"`` or ``"null"``.
                - ``theta_deg`` (float): Target elevation. Units: degrees.
                - ``phi_deg`` (float): Target azimuth. Units: degrees.
                - ``gain_dbi`` (float): Directivity at the nearest grid point. Units: dBi.
                - ``null_depth_db`` (float | None): gain_dbi − global_peak_dbi (≤ 0).
                  Present only for null directives; None for peak directives.
                - ``cost_term`` (float): lambda_k × C_k, signed as in the cost function.

            - ``peak_to_null_ratio_db`` (float | None): Difference in gain_dbi between
              the first peak and first null directive. None if either type is absent.
    """
    # ── Array factor ──────────────────────────────────────────────
    # Coherent superposition: AF(theta, phi) = sum_n w_n * E_n(theta, phi)
    array_factor = compute_array_factor(weights_complex, element_patterns_stacked)
    # [MATLAB] array_factor = sum(permute(weights_complex, [1 3 4 2]) .* element_patterns_stacked, 1);

    # Power density over the full angular grid.
    power_linear = np.abs(array_factor) ** 2
    # [MATLAB] power_linear = abs(array_factor) .^ 2;

    # Directivity grid: 4π · |AF|² / ∫|AF|² dΩ — absolute V/m scale cancels.
    directivity_dbi_grid = _compute_directivity_dbi_grid(array_factor, theta_deg, phi_deg)
    global_peak_dbi      = float(np.max(directivity_dbi_grid))
    # [MATLAB] global_peak_dbi = max(directivity_dbi_grid(:));

    # ── Per-directive evaluation ───────────────────────────────────
    directive_metrics  = []
    total_cost         = 0.0
    first_peak_dbi     = None
    first_null_dbi     = None

    for directive in directives:
        directive_type   = directive["type"]
        target_theta_deg = directive["theta"]
        target_phi_deg   = directive.get("phi",    DEFAULT_TARGET_PHI_DEG)
        width_deg        = directive["width"]
        directive_weight = directive.get("weight", DEFAULT_DIRECTIVE_WEIGHT)

        # Directivity at the nearest grid point to the target angle.
        theta_idx = _nearest_index(theta_deg, target_theta_deg)
        phi_idx   = _nearest_index(phi_deg,   target_phi_deg)
        gain_dbi  = float(directivity_dbi_grid[theta_idx, phi_idx])
        # [MATLAB] [~, ti] = min(abs(theta_deg - target_theta_deg));
        # [MATLAB] [~, pi] = min(abs(phi_deg   - target_phi_deg));
        # [MATLAB] gain_dbi = directivity_dbi_grid(ti, pi);

        # Mean power in the angular window — used for the cost term (matches cost function).
        mask              = angular_window_mask(
            theta_deg, phi_deg, target_theta_deg, target_phi_deg, width_deg
        )
        mean_window_power = np.mean(power_linear[mask])
        # [MATLAB] mean_window_power = mean(power_linear(mask));

        if directive_type == "peak":
            cost_term     = -directive_weight * float(mean_window_power)
            null_depth_db = None
            if first_peak_dbi is None:
                first_peak_dbi = gain_dbi
        else:  # "null"
            cost_term     = +directive_weight * float(mean_window_power)
            # Null depth: how far below the global peak this null sits.
            null_depth_db = gain_dbi - global_peak_dbi
            # [MATLAB] null_depth_db = gain_dbi - global_peak_dbi;
            if first_null_dbi is None:
                first_null_dbi = gain_dbi

        total_cost += cost_term

        directive_metrics.append({
            "type":          directive_type,
            "theta_deg":     float(target_theta_deg),
            "phi_deg":       float(target_phi_deg),
            "gain_dbi":      gain_dbi,
            "null_depth_db": float(null_depth_db) if null_depth_db is not None else None,
            "cost_term":     float(cost_term),
        })

    # ── Peak-to-null ratio ─────────────────────────────────────────
    if first_peak_dbi is not None and first_null_dbi is not None:
        peak_to_null_ratio_db = first_peak_dbi - first_null_dbi
    else:
        peak_to_null_ratio_db = None

    return {
        "total_cost":            float(total_cost),
        "iteration_count":       len(cost_history),
        "global_peak_dbi":       global_peak_dbi,
        "directive_metrics":     directive_metrics,
        "peak_to_null_ratio_db": (float(peak_to_null_ratio_db)
                                  if peak_to_null_ratio_db is not None else None),
    }
