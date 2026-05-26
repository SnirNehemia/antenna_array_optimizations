# ══════════════════════════════════════════════════════════════════
# COST_FUNCTION
# Composite cost function J(x) for antenna array pattern optimization.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────
import numpy as np

# ────────────────────────── CONSTANTS ─────────────────────────────

# Default directive weight when not specified by the user.
DEFAULT_DIRECTIVE_WEIGHT = 1.0

# Default phi target when not specified by a directive.
DEFAULT_TARGET_PHI_DEG = 0.0

# Valid optimization modes passed to build_cost_function.
VALID_MODES = ("standard", "phase_only", "amplitude_only")

# Valid directive types.
VALID_DIRECTIVE_TYPES = ("peak", "null")


# ────────────────────────── VARIABLE ENCODING ─────────────────────

def x_to_weights(x):
    """Decode the 2N real optimization variable vector into N complex weights.

    The encoding interleaves real and imaginary parts:
    x = [Re(w_0), Im(w_0), Re(w_1), Im(w_1), ..., Re(w_{N-1}), Im(w_{N-1})]

    Args:
        x (np.ndarray): Real optimization variable vector, shape (2N,).
            Units: dimensionless (V/V).

    Returns:
        np.ndarray: Complex weight vector, shape (N,). Units: dimensionless (V/V).
    """
    weights_real = x[0::2]  # [MATLAB] x(1:2:end)
    weights_imag = x[1::2]  # [MATLAB] x(2:2:end)
    return weights_real + 1j * weights_imag


def weights_to_x(weights_complex):
    """Encode N complex weights into the 2N real optimization variable vector.

    Interleaves real and imaginary parts so that element n occupies positions
    2n (real) and 2n+1 (imaginary) in the returned vector.

    Args:
        weights_complex (np.ndarray): Complex weight vector, shape (N,).
            Units: dimensionless (V/V).

    Returns:
        np.ndarray: Real optimization variable vector, shape (2N,).
            Units: dimensionless (V/V).
    """
    # Stack real and imaginary columns side-by-side, then flatten row-major.
    x = np.column_stack([weights_complex.real, weights_complex.imag]).ravel()
    # [MATLAB] x = reshape([real(weights_complex), imag(weights_complex)]', 1, []);
    return x


# ────────────────────────── ARRAY FACTOR ──────────────────────────

def compute_array_factor(weights_complex, element_patterns_stacked):
    """Compute the coherent array factor over the full angular grid.

    Implements AF(theta, phi) = sum_{n=0}^{N-1} w_n * E_n(theta, phi),
    the coherent superposition of all weighted element patterns.

    Args:
        weights_complex (np.ndarray): Complex element weights, shape (N_elements,).
            Units: dimensionless (V/V).
        element_patterns_stacked (np.ndarray): Complex element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m.

    Returns:
        np.ndarray: Complex array factor, shape (N_theta, N_phi). Units: V/m.
    """
    # Broadcast weights across the spatial grid dimensions before summing.
    # w_n * E_n broadcasts: (N,1,1) * (N, N_theta, N_phi) → (N, N_theta, N_phi)
    array_factor = np.sum(
        weights_complex[:, np.newaxis, np.newaxis] * element_patterns_stacked,
        axis=0,
    )
    # [MATLAB] for n = 1:N_elements
    # [MATLAB]     array_factor = array_factor + weights_complex(n) * squeeze(element_patterns_stacked(n,:,:));
    # [MATLAB] end
    return array_factor


# ────────────────────────── ANGULAR WINDOW ────────────────────────

def angular_window_mask(theta_deg, phi_deg,
                        target_theta_deg, target_phi_deg,
                        theta_width_deg, phi_width_deg):
    """Build a boolean 2D mask covering an angular window on the (theta, phi) grid.

    The window is a rectangular region centred on (target_theta, target_phi).
    Elevation and azimuth half-widths are specified independently so the window
    can be asymmetric (e.g. narrow in phi, wide in theta).

    Args:
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        target_theta_deg (float): Window centre in elevation. Units: degrees.
        target_phi_deg (float): Window centre in azimuth. Units: degrees.
        theta_width_deg (float): Full angular width in elevation. Units: degrees.
        phi_width_deg (float): Full angular width in azimuth. Units: degrees.

    Returns:
        np.ndarray: Boolean mask, shape (N_theta, N_phi). True where the grid
            point falls inside the angular window.
    """
    theta_hw = theta_width_deg / 2.0
    phi_hw   = phi_width_deg   / 2.0

    # Identify which theta and phi values fall within the window.
    theta_in_window = np.abs(theta_deg - target_theta_deg) <= theta_hw
    # [MATLAB] theta_in_window = abs(theta_deg - target_theta_deg) <= theta_hw;
    phi_in_window = np.abs(phi_deg - target_phi_deg) <= phi_hw
    # [MATLAB] phi_in_window = abs(phi_deg - target_phi_deg) <= phi_hw;

    # Outer product: mask[i, j] = True iff theta_i and phi_j are both in window.
    mask = theta_in_window[:, np.newaxis] & phi_in_window[np.newaxis, :]
    # [MATLAB] mask = theta_in_window' & phi_in_window;
    return mask


# ────────────────────────── DIRECTIVE COST ────────────────────────

def _directive_cost(power_grid, mask, directive_type, sin_theta):
    """Compute the scalar cost contribution for a single beam-shaping directive.

    All cost terms are formulated as quantities to **minimize** so that
    L-BFGS-B (a minimizer) drives the array pattern toward the desired shape.

    - Peak:  C = -mean_sa(|AF|²) in window  (negative → minimizing increases gain)
    - Null:  C = +mean_sa(|AF|²) in window  (positive → minimizing suppresses gain)

    The window mean is weighted by solid angle (sin θ · Δθ · Δφ) so that pixels
    near the poles (where sin θ ≈ 0) do not over-contribute relative to their
    actual angular area.

    Args:
        power_grid (np.ndarray): |AF|² over the full grid, shape (N_theta, N_phi).
            Units: (V/m)². Pre-computed by the caller once per cost evaluation.
        mask (np.ndarray): Boolean angular window mask, shape (N_theta, N_phi).
        directive_type (str): Either ``"peak"`` or ``"null"``.
        sin_theta (np.ndarray): sin(θ) for each elevation sample, shape (N_theta,).
            Used as the solid-angle weighting factor.

    Returns:
        float: Scalar cost contribution for this directive. Units: (V/m)².

    Raises:
        ValueError: If ``directive_type`` is not ``"peak"`` or ``"null"``.
    """
    if directive_type not in VALID_DIRECTIVE_TYPES:
        raise ValueError(
            f"Unknown directive type '{directive_type}'. "
            f"Valid types: {VALID_DIRECTIVE_TYPES}."
        )

    # Solid-angle-weighted mean power inside the window.
    # weights_sa[i,j] = sin(θ_i) when mask[i,j] is True, 0 otherwise.
    # [MATLAB] weights_sa = mask .* sin(deg2rad(theta_deg))';
    weights_sa   = mask * sin_theta[:, np.newaxis]
    total_weight = float(np.sum(weights_sa))
    mean_power   = float(np.sum(power_grid * weights_sa)) / max(total_weight, 1e-30)
    # [MATLAB] mean_power = sum(power_grid(:) .* weights_sa(:)) / max(sum(weights_sa(:)), 1e-30);

    if directive_type == "peak":
        # Maximizing gain = minimizing the negative of mean power in the beam window.
        return -mean_power
    else:  # "null"
        # Suppressing gain = minimizing mean power in the null window.
        return mean_power


# ────────────────────────── COST FUNCTION BUILDER ─────────────────

def build_cost_function(element_patterns_stacked, theta_deg, phi_deg,
                        directives, mode="standard"):
    """Build and return the composite cost function J(x) as a callable.

    Pre-computes angular window masks for all directives (constant throughout
    optimization) and returns a closure that evaluates J(x) on every call.

    The composite cost is:
        J(x) = sum_k  lambda_k * C_k(x)
    where lambda_k is the user-supplied ``weight`` for directive k and C_k is the
    per-directive cost (peak or null term).

    Args:
        element_patterns_stacked (np.ndarray): Complex element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        directives (list[dict]): Beam-shaping directives. Each dict must have:
            - ``type`` (str): ``"peak"`` or ``"null"``.
            - ``theta`` (float): Target elevation angle. Units: degrees.
            - ``width`` (float): Symmetric angular window width (used when
              ``theta_width`` / ``phi_width`` are absent). Units: degrees.
            Optional keys:
            - ``theta_width`` (float): Elevation window width (overrides ``width``).
            - ``phi_width`` (float): Azimuth window width (overrides ``width``).
            - ``phi`` (float): Target azimuth angle (default 0.0). Units: degrees.
            - ``weight`` (float): Directive weight lambda_k (default 1.0).
        mode (str): Optimization mode. One of:

            - ``"standard"``: x has length 2N; weights decoded as complex Re/Im pairs.
            - ``"phase_only"``: x has length 2N; weights are normalized to unit amplitude.
            - ``"amplitude_only"``: x has length N; weights are real (phase = 0).

    Returns:
        callable: A function ``cost_fn(x) -> float`` that evaluates J(x).

    Raises:
        ValueError: If ``mode`` is not one of the valid modes, or if any directive
            has an unknown ``type``.
    """
    if mode not in VALID_MODES:
        raise ValueError(
            f"Unknown optimization mode '{mode}'. Valid modes: {VALID_MODES}."
        )

    # Pre-compute angular window masks — these are fixed for the entire optimization.
    # theta_width / phi_width are independent; "width" is the symmetric fallback.
    directive_masks = []
    for directive in directives:
        target_theta_deg = directive["theta"]
        target_phi_deg   = directive.get("phi", DEFAULT_TARGET_PHI_DEG)
        sym_width        = directive.get("width", None)
        theta_width_deg  = directive.get("theta_width", sym_width)
        phi_width_deg    = directive.get("phi_width",   sym_width)
        mask = angular_window_mask(
            theta_deg, phi_deg, target_theta_deg, target_phi_deg,
            theta_width_deg, phi_width_deg,
        )
        directive_masks.append(mask)

    # Pre-compute solid-angle weights: sin(θ) for each elevation sample.
    # Captured by the cost_fn closure; computed once per build_cost_function call.
    # [MATLAB] sin_theta = sin(deg2rad(theta_deg));
    sin_theta = np.sin(np.deg2rad(theta_deg))

    def cost_fn(x):
        """Evaluate the composite cost J(x) for a given variable vector x."""
        # Decode x into complex weights according to the optimization mode.
        if mode == "amplitude_only":
            # Amplitude-only: x is a real amplitude vector of length N.
            # Phase is fixed at 0 for all elements.
            weights_complex = x.astype(complex)
            # [MATLAB] weights_complex = complex(x, zeros(size(x)));
        else:
            # Standard and phase-only: x is a 2N real vector (Re/Im interleaved).
            weights_complex = x_to_weights(x)

            if mode == "phase_only":
                # Enforce unit amplitude: project each weight onto the unit circle.
                # np.where avoids division by zero for any zero-initialized weights.
                amplitudes = np.abs(weights_complex)
                safe_amplitudes = np.where(amplitudes > 0.0, amplitudes, 1.0)
                # [MATLAB] safe_amplitudes = max(abs(weights_complex), eps);
                weights_complex = weights_complex / safe_amplitudes

        # Power-normalize weights so the cost is independent of global amplitude scale.
        # Physical model: fixed total source power split across elements → fair comparison
        # between different amplitude distributions (power-splitter constraint).
        # AF_norm = AF / ||w||₂ ; cost ∝ |AF_norm|² = |AF|² / Σ|w_n|²
        # [MATLAB] w_power = sum(abs(weights_complex).^2); weights_complex = weights_complex / sqrt(max(w_power, 1e-30));
        w_power    = float(np.sum(np.abs(weights_complex) ** 2))
        safe_power = max(w_power, 1e-30)
        weights_complex = weights_complex / np.sqrt(safe_power)

        # Array factor: coherent superposition of power-normalized element patterns.
        array_factor = compute_array_factor(weights_complex, element_patterns_stacked)

        # Pre-compute power grid once; reused for every directive this call.
        # [MATLAB] power_grid = abs(array_factor) .^ 2;
        power_grid = np.abs(array_factor) ** 2

        # Composite cost: weighted sum over all directives.
        total_cost = 0.0
        for directive, mask in zip(directives, directive_masks):
            directive_type   = directive["type"]
            directive_weight = directive.get("weight", DEFAULT_DIRECTIVE_WEIGHT)
            term = _directive_cost(power_grid, mask, directive_type, sin_theta)
            total_cost += directive_weight * term

        return float(total_cost)

    return cost_fn
