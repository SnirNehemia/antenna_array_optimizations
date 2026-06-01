# ══════════════════════════════════════════════════════════════════
# METRICS
# Post-optimization evaluation: per-directive gain scores and cost breakdown.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────
import numpy as np

from src.cost.cost_function import compute_array_factor, build_directive_physical_masks

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

def _compute_directivity_dbi_grid(array_factor, theta_deg, phi_deg,
                                   normalizer_power=None):
    """Compute directivity in dBi over the full (N_theta × N_phi) angular grid.

    Uses the discrete rectangle-rule approximation to the spherical integral:

        P_total ≈ Σ_θ Σ_φ |AF(θ,φ)|² · sin(θ) · Δθ · Δφ   [(V/m)² · sr]
        D(θ,φ)  =  4π · |AF(θ,φ)|² / P_total               [dimensionless]

    The absolute V/m² scale of the element patterns cancels exactly in the
    ratio, so the result is independent of CST normalisation.

    Following the CST / IEEE Std 149 "partial directivity" convention, callers
    that have both polarisation stacks should pass the combined radiated power
    ``P_copol + P_cross`` via ``normalizer_power`` so that copol and cross
    directivities share the same denominator and their sum equals the total
    directivity at every point.

    Args:
        array_factor (np.ndarray): Complex array factor, shape (N_theta, N_phi).
            Units: V/m (arbitrary absolute scale).
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
            Must cover 0° – 180° for a full-sphere integral.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        normalizer_power (float | None): Pre-computed total radiated power used as
            the denominator instead of the self-computed power of ``array_factor``.
            Pass ``P_copol + P_cross`` here when both polarisation stacks are
            available to get CST-convention partial directivity. Default: None
            (use the power radiated by ``array_factor`` alone).

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

    if normalizer_power is not None:
        p_total = float(normalizer_power)
    else:
        # Numerical spherical integral: Σ |AF|² sin(θ) Δθ Δφ
        p_total = np.sum(power_grid * sin_theta[:, np.newaxis]) * dtheta_rad * dphi_rad
    # [MATLAB] p_total = sum(sum(power_grid .* sin(theta_rad)')) * dtheta_rad * dphi_rad;

    directivity_linear = (4.0 * np.pi * power_grid) / max(p_total, 1e-30)
    # [MATLAB] directivity_linear = 4 * pi * power_grid / p_total;

    return 10.0 * np.log10(directivity_linear + LOG10_EPSILON)
    # [MATLAB] directivity_dbi_grid = 10 * log10(directivity_linear + eps);


# ────────────────────────── HPBW HELPER ───────────────────────────

def _compute_hpbw(cut_dbi, peak_idx, opposite_cut=None):
    """3 dB beamwidth by scanning outward from peak_idx with linear interpolation.

    When the scan reaches a grid boundary without finding a crossing and
    ``opposite_cut`` is provided (the θ-cut at φ+180°), the beam is assumed to
    wrap around the pole and the search continues in that cut.  This corrects
    the HPBW for beams pointing near θ=0° or θ=180°.

    Returns the width in grid-index units; the caller multiplies by the angular
    step size to convert to degrees.

    Args:
        cut_dbi (np.ndarray): 1-D directivity cut in dBi.
        peak_idx (int): Index of the global peak within ``cut_dbi``.
        opposite_cut (np.ndarray | None): θ-cut at azimuth φ+180° for pole
            wrap-around correction.  None disables wrap-around (φ-cut usage).

    Returns:
        float: Full 3 dB beamwidth in grid-index units.
    """
    threshold = float(cut_dbi[peak_idx]) - 3.0
    n = len(cut_dbi)

    # ── Left (decreasing-index) crossing ──────────────────────────
    # cut[i] < threshold ≤ cut[i+1] at crossing; denom > 0, t ∈ [0, 1].
    left = float(peak_idx)
    left_found = False
    for i in range(peak_idx - 1, -1, -1):
        if cut_dbi[i] < threshold:
            denom = float(cut_dbi[i + 1] - cut_dbi[i])   # > 0
            t = (threshold - float(cut_dbi[i])) / (denom if abs(denom) > 1e-30 else 1e-30)
            left = float(i) + float(np.clip(t, 0.0, 1.0))
            left_found = True
            break

    if not left_found and opposite_cut is not None:
        # Beam extends past θ=0° (north pole).  Continue in the opposite-azimuth
        # cut, scanning outward from index 0.  The virtual left crossing is
        # opp_extra steps before the pole (index −opp_extra in extended coords).
        for i in range(1, n):
            if opposite_cut[i] < threshold:
                denom = float(opposite_cut[i] - opposite_cut[i - 1])   # < 0
                t = (threshold - float(opposite_cut[i - 1])) / (denom if abs(denom) > 1e-30 else -1e-30)
                opp_extra = float(i - 1) + float(np.clip(t, 0.0, 1.0))
                left = -opp_extra
                break

    # ── Right (increasing-index) crossing ─────────────────────────
    # cut[i-1] ≥ threshold > cut[i] at crossing; denom < 0; must use signed.
    right = float(peak_idx)
    right_found = False
    for i in range(peak_idx + 1, n):
        if cut_dbi[i] < threshold:
            denom = float(cut_dbi[i] - cut_dbi[i - 1])   # < 0
            t = (threshold - float(cut_dbi[i - 1])) / (denom if abs(denom) > 1e-30 else -1e-30)
            right = float(i - 1) + float(np.clip(t, 0.0, 1.0))
            right_found = True
            break

    if not right_found and opposite_cut is not None:
        # Beam extends past θ=180° (south pole).  Continue in the opposite-azimuth
        # cut, scanning backward from index n-1.  The virtual right crossing is
        # opp_extra steps past the pole (index (n-1) + opp_extra in extended coords).
        for i in range(n - 2, -1, -1):
            if opposite_cut[i] < threshold:
                denom = float(opposite_cut[i + 1] - opposite_cut[i])   # > 0
                t = (threshold - float(opposite_cut[i])) / (denom if abs(denom) > 1e-30 else 1e-30)
                opp_extra = float(n - 1) - (float(i) + float(np.clip(t, 0.0, 1.0)))
                right = float(n - 1) + opp_extra
                break

    return right - left


# ────────────────────────── PUBLIC INTERFACE ──────────────────────

def evaluate_metrics(element_patterns_stacked, theta_deg, phi_deg,
                     weights_complex, directives, cost_history,
                     precomputed_array_factor=None,
                     normalizer_power=None):
    """Evaluate post-optimization performance metrics for all beam-shaping directives.

    Computes the array factor from the final optimized weights, then for each
    directive reports the realized power level at the target angle (in dB) and
    the cost contribution. Null depth is expressed relative to the global pattern peak.

    Args:
        element_patterns_stacked (np.ndarray): Complex element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m. Ignored when
            ``precomputed_array_factor`` is provided.
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
        precomputed_array_factor (np.ndarray | None): Optional pre-built array
            factor, shape (N_theta, N_phi). When provided, ``compute_array_factor``
            is skipped — used by the ``total`` polarization mode in manual_weights.py
            where the displayed power is the orthogonal sum
            ``|AF_copol|² + |AF_xpol|²`` and cannot be expressed as a single
            coherent superposition of one stacked pattern.
        normalizer_power (float | None): Total radiated power used as the
            directivity denominator (passed through to
            ``_compute_directivity_dbi_grid``). When both polarisation stacks are
            available pass ``P_copol + P_cross`` here for CST-convention partial
            directivity (D_copol + D_cross = D_total). Default: None (use the
            power radiated by the array factor alone).

    Returns:
        dict: Metrics summary with keys:

            - ``total_cost`` (float): Composite cost J(x*) at the optimized weights.
            - ``iteration_count`` (int): Number of optimizer iterations recorded.
            - ``global_peak_dbi`` (float): Peak directivity over the full grid. Units: dBi.
            - ``directive_metrics`` (list[dict]): Per-directive results. Each entry has:

                - ``type`` (str): ``"peak"`` or ``"null"``.
                - ``theta_deg`` (float): Target elevation. Units: degrees.
                - ``phi_deg`` (float): Target azimuth. Units: degrees.
                - ``gain_dbi`` (float): Peak directivity inside the directive's
                  angular window — the strongest point of the on-screen window box.
                  Resolves pole-crossing (θ<0° or θ>180°) and φ wrap-around. Units: dBi.
                - ``null_depth_db`` (float | None): gain_dbi − global_peak_dbi (≤ 0).
                  Present only for null directives; None for peak directives.
                - ``cost_term`` (float): lambda_k × C_k, signed as in the cost function.

            - ``peak_to_null_ratio_db`` (float | None): Difference in gain_dbi between
              the first peak and first null directive. None if either type is absent.
    """
    # ── Array factor ──────────────────────────────────────────────
    # Coherent superposition: AF(theta, phi) = sum_n w_n * E_n(theta, phi).
    # If the caller has already computed (or synthesised) the array factor —
    # e.g. for the total-polarization power-sum case — use it directly.
    if precomputed_array_factor is not None:
        array_factor = precomputed_array_factor
    else:
        array_factor = compute_array_factor(weights_complex, element_patterns_stacked)
    # [MATLAB] array_factor = sum(permute(weights_complex, [1 3 4 2]) .* element_patterns_stacked, 1);

    # Power density over the full angular grid.
    power_linear = np.abs(array_factor) ** 2
    # [MATLAB] power_linear = abs(array_factor) .^ 2;

    # Directivity grid: 4π · |AF|² / P_total — absolute V/m scale cancels.
    # normalizer_power, when provided, is the combined P_copol + P_cross so that
    # copol and cross partial directivities share the same denominator (CST convention).
    directivity_dbi_grid = _compute_directivity_dbi_grid(
        array_factor, theta_deg, phi_deg, normalizer_power=normalizer_power
    )
    global_peak_dbi      = float(np.max(directivity_dbi_grid))
    # [MATLAB] global_peak_dbi = max(directivity_dbi_grid(:));

    # ── Global peak location and 3 dB HPBW ───────────────────────
    theta_peak_idx, phi_peak_idx = np.unravel_index(
        np.argmax(directivity_dbi_grid), directivity_dbi_grid.shape
    )
    global_peak_theta_deg = float(theta_deg[theta_peak_idx])
    global_peak_phi_deg   = float(phi_deg[phi_peak_idx])

    # θ-HPBW: 1-D cut at peak φ, with pole wrap-around via opposite azimuth.
    # When the beam points near θ=0° or θ=180°, one side of the 3 dB contour
    # extends past the pole; the opposite-azimuth cut (φ+180°) continues it.
    n_phi         = len(phi_deg)
    theta_step    = float(theta_deg[1] - theta_deg[0]) if len(theta_deg) > 1 else 1.0
    theta_cut     = directivity_dbi_grid[:, phi_peak_idx]
    opp_phi_idx   = (int(phi_peak_idx) + n_phi // 2) % n_phi
    opp_theta_cut = directivity_dbi_grid[:, opp_phi_idx]
    hpbw_theta_deg = _compute_hpbw(theta_cut, int(theta_peak_idx),
                                   opposite_cut=opp_theta_cut) * theta_step

    # φ-HPBW: 1-D cut at peak θ, rolled so peak is at centre (handles wrap-around)
    phi_step    = float(phi_deg[1] - phi_deg[0]) if n_phi > 1 else 1.0
    phi_cut_raw = directivity_dbi_grid[theta_peak_idx, :]
    roll_shift  = n_phi // 2 - int(phi_peak_idx)
    phi_cut     = np.roll(phi_cut_raw, roll_shift)
    hpbw_phi_deg = _compute_hpbw(phi_cut, n_phi // 2) * phi_step

    # ── Per-directive evaluation ───────────────────────────────────
    # Pre-compute solid-angle weights: sin(θ) broadcast to (N_theta, 1).
    # Reused for every directive's window mean to match the cost function weighting.
    # [MATLAB] sin_theta_grid = sin(deg2rad(theta_deg))';  % column vector, broadcast over phi
    sin_theta_grid = np.sin(np.deg2rad(theta_deg))[:, np.newaxis]

    # Physical-grid window masks — identical to the cost function and the GUI
    # overlay box. These resolve pole-crossing (θ<0° or θ>180°) and φ 0°/360°
    # wrap-around, so a directive at e.g. θ=−30° is mapped to its physical mirror
    # location (θ=30°, φ+180°) instead of being clamped to the grid edge.
    # [MATLAB] phys_masks = build_directive_physical_masks(theta_deg, phi_deg, directives);
    phys_masks = build_directive_physical_masks(theta_deg, phi_deg, directives)

    directive_metrics  = []
    total_cost         = 0.0
    first_peak_dbi     = None
    first_null_dbi     = None

    for directive, mask in zip(directives, phys_masks):
        directive_type   = directive["type"]
        target_theta_deg = directive["theta"]
        target_phi_deg   = directive.get("phi",    DEFAULT_TARGET_PHI_DEG)
        directive_weight = directive.get("weight", DEFAULT_DIRECTIVE_WEIGHT)

        # Realized directivity: the strongest point inside the directive's physical
        # window. For a peak this is the achieved gain; for a null it is the
        # worst-case leakage — both match the brightest pixel of the on-screen box.
        # Falls back to the mapped-centre nearest grid point only if the window
        # somehow covers no grid points (degenerate width).
        if mask.any():
            gain_dbi = float(np.max(directivity_dbi_grid[mask]))
            # [MATLAB] gain_dbi = max(directivity_dbi_grid(mask));
        else:
            theta_idx = _nearest_index(theta_deg, target_theta_deg)
            phi_idx   = _nearest_index(phi_deg,   target_phi_deg)
            gain_dbi  = float(directivity_dbi_grid[theta_idx, phi_idx])

        # Solid-angle-weighted mean power over the same physical window — matches
        # the cost function's default "mean" aggregation.
        # Near-pole pixels (sin θ ≈ 0) are de-emphasized proportional to their solid angle.
        # [MATLAB] w = mask .* sin_theta_grid; mean_window_power = sum(power_linear(:).*w(:)) / max(sum(w(:)), 1e-30);
        weights_sa        = mask.astype(float) * sin_theta_grid
        total_weight      = float(np.sum(weights_sa))
        mean_window_power = (float(np.sum(power_linear * weights_sa)) / max(total_weight, 1e-30)
                             if total_weight > 1e-30 else 0.0)

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
        "global_peak_theta_deg": global_peak_theta_deg,
        "global_peak_phi_deg":   global_peak_phi_deg,
        "hpbw_theta_deg":        hpbw_theta_deg,
        "hpbw_phi_deg":          hpbw_phi_deg,
        "directive_metrics":     directive_metrics,
        "peak_to_null_ratio_db": (float(peak_to_null_ratio_db)
                                  if peak_to_null_ratio_db is not None else None),
    }
