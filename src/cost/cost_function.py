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

# Valid cost aggregation methods per directive.
VALID_AGGREGATIONS = ("mean", "max", "min")


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

    This function operates on whatever theta/phi arrays are passed in.  Callers
    that need pole-crossing or 0°/360° wrap-around should pass the pre-built
    extended grids from build_cost_function rather than the raw physical grids.

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

def _directive_cost(masked_power, sin_weights, directive_type, aggregation="mean"):
    """Compute the scalar cost contribution for a single beam-shaping directive.

    Operates on the sparse set of angular-window sample points (pre-extracted
    by build_cost_function from the extended grid).  All cost terms are
    formulated as quantities to **minimize**.

    - Peak:  C = -metric  (negative → minimizing increases gain)
    - Null:  C = +metric  (positive → minimizing suppresses gain)

    Args:
        masked_power (np.ndarray): |AF|² at each in-window sample, shape (K,).
            Units: (V/m)².
        sin_weights (np.ndarray): |sin(θ)| solid-angle weight per sample, shape (K,).
            Used only when aggregation == "mean".
        directive_type (str): ``"peak"`` or ``"null"``.
        aggregation (str): ``"mean"``, ``"max"``, or ``"min"``.

    Returns:
        float: Scalar cost contribution. Units: (V/m)².

    Raises:
        ValueError: If ``directive_type`` or ``aggregation`` is invalid.
    """
    if directive_type not in VALID_DIRECTIVE_TYPES:
        raise ValueError(
            f"Unknown directive type '{directive_type}'. "
            f"Valid types: {VALID_DIRECTIVE_TYPES}."
        )
    if aggregation not in VALID_AGGREGATIONS:
        raise ValueError(
            f"Unknown aggregation '{aggregation}'. "
            f"Valid options: {VALID_AGGREGATIONS}."
        )

    if len(masked_power) == 0:
        return 0.0

    if aggregation == "mean":
        # Solid-angle-weighted mean power.
        total_w = float(np.sum(sin_weights))
        metric  = float(np.dot(masked_power, sin_weights)) / max(total_w, 1e-30)
        # [MATLAB] metric = sum(masked_power .* sin_weights) / max(sum(sin_weights), 1e-30);
    elif aggregation == "max":
        # Worst-case (maximum) power in window.
        # Best for nulls: drives down the highest sidelobe inside the window.
        metric = float(np.max(masked_power))
        # [MATLAB] metric = max(masked_power);
    else:  # "min"
        # Best-case (minimum) power in window.
        # Best for peaks: enforces a flat beam by maximising the lowest point.
        metric = float(np.min(masked_power))
        # [MATLAB] metric = min(masked_power);

    if directive_type == "peak":
        return -metric
    else:  # "null"
        return +metric



# ────────────────────────── PHYSICAL MASK BUILDER ─────────────────

def build_directive_physical_masks(theta_deg, phi_deg, directives):
    """Build the physical-grid boolean mask for each directive.

    Uses the same extended-grid logic as build_cost_function so the returned
    masks exactly match the angular windows used during optimization, including
    phi 0°/360° wrap-around and theta pole-crossing at both poles.

    Args:
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        directives (list[dict]): Beam-shaping directives (same format as
            build_cost_function).

    Returns:
        list[np.ndarray]: One boolean (N_theta, N_phi) array per directive.
            True where the directive's angular window covers that physical grid point.
    """
    N_theta = len(theta_deg)
    N_phi   = len(phi_deg)
    PHI_180_SHIFT = N_phi // 2

    theta_ext_deg = np.concatenate([
        -theta_deg[1:][::-1],
        theta_deg,
        360.0 - theta_deg[N_theta - 2::-1],
    ])
    phi_ext_deg = np.concatenate([phi_deg - 360.0, phi_deg, phi_deg + 360.0])

    theta_phys_idx = np.concatenate([
        np.arange(N_theta - 1, 0, -1),
        np.arange(N_theta),
        np.arange(N_theta - 2, -1, -1),
    ])
    phi_offset = np.concatenate([
        np.full(N_theta - 1, PHI_180_SHIFT, dtype=int),
        np.zeros(N_theta, dtype=int),
        np.full(N_theta - 1, PHI_180_SHIFT, dtype=int),
    ])

    masks = []
    for directive in directives:
        target_theta_deg = directive["theta"]
        target_phi_deg   = directive.get("phi", DEFAULT_TARGET_PHI_DEG)
        sym_width        = directive.get("width", None)
        theta_width_deg  = directive.get("theta_width", sym_width)
        phi_width_deg    = directive.get("phi_width",   sym_width)

        ext_mask = angular_window_mask(
            theta_ext_deg, phi_ext_deg,
            target_theta_deg, target_phi_deg,
            theta_width_deg, phi_width_deg,
        )
        ext_t, ext_p = np.nonzero(ext_mask)
        phys_theta = theta_phys_idx[ext_t]
        phys_phi   = (ext_p % N_phi + phi_offset[ext_t]) % N_phi
        # [MATLAB] phys_mask = false(N_theta, N_phi); phys_mask(sub2ind([N_theta,N_phi], phys_theta+1, phys_phi+1)) = true;
        phys_mask = np.zeros((N_theta, N_phi), dtype=bool)
        phys_mask[phys_theta, phys_phi] = True
        masks.append(phys_mask)

    return masks


# ────────────────────────── COST FUNCTION BUILDER ─────────────────

def build_cost_function(element_patterns_stacked, theta_deg, phi_deg,
                        directives, mode="standard",
                        element_patterns_secondary=None):
    """Build and return the composite cost function J(x) as a callable.

    Constructs an extended (θ, φ) grid that mirrors the physical pattern at the
    north pole (θ<0°) and south pole (θ>180°) and tiles it three times in φ.
    Angular window masks are computed on this extended grid so that rectangular
    windows near the poles or straddling the φ=0°/360° seam are handled
    correctly without special-case logic.

    Each directive's mask is then reduced to a sparse list of physical (θ, φ)
    indices so the closure only touches the relevant grid points on every call.

    The composite cost is:
        J(x) = sum_k  lambda_k * C_k(x)
    where lambda_k is the user-supplied ``weight`` for directive k and C_k is the
    per-directive cost (peak or null term).

    Args:
        element_patterns_stacked (np.ndarray): Complex element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
            Expected range: 0° to 180°, uniformly spaced.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
            Expected range: 0° to <360°, uniformly spaced, N_phi even.
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
            - ``aggregation`` (str): ``"mean"`` | ``"max"`` | ``"min"`` (default ``"mean"``).
        mode (str): Optimization mode. One of:

            - ``"standard"``: x has length 2N; weights decoded as complex Re/Im pairs.
            - ``"phase_only"``: x has length 2N; weights are normalized to unit amplitude.
            - ``"amplitude_only"``: x has length N; weights are real (phase = 0).
        element_patterns_secondary (np.ndarray | None): Optional second polarisation
            stack, shape (N_elements, N_theta, N_phi). When provided, the cost
            function operates on total field power:
            ``|AF_primary|² + |AF_secondary|²``.
            Pass ``cross_complex`` here for ``polarization: "total"`` mode.
            Default: None (single-stack, single-polarisation).

    Returns:
        callable: A function ``cost_fn(x) -> float`` that evaluates J(x).

    Raises:
        ValueError: If ``mode`` is not one of the valid modes, or if any directive
            has an unknown ``type`` or ``aggregation``.
    """
    if mode not in VALID_MODES:
        raise ValueError(
            f"Unknown optimization mode '{mode}'. Valid modes: {VALID_MODES}."
        )

    N_theta = len(theta_deg)
    N_phi   = len(phi_deg)
    # Assumes a uniform phi grid with N_phi points covering 360° (step = 360°/N_phi).
    # PHI_180_SHIFT is the number of phi columns that corresponds to a 180° rotation.
    PHI_180_SHIFT = N_phi // 2  # [MATLAB] PHI_180_SHIFT = N_phi / 2;

    # ── Extended grids ────────────────────────────────────────────────
    # Theta is mirrored at both poles so windows near θ=0° or θ=180° wrap
    # correctly into the back hemisphere.
    #
    # theta_ext layout (3N-2 points):
    #   top mirror  : -180°, ..., -1°  → physical (|θ|, φ+180°)
    #   main        :    0°, ..., 180° → physical (θ, φ)
    #   bottom mirror: 181°, ..., 360° → physical (360°-θ, φ+180°)
    #
    # phi_ext layout (3*N_phi points):
    #   left tile : φ-360°  → physical φ (mod N_phi)
    #   centre    : φ       → physical φ
    #   right tile: φ+360°  → physical φ (mod N_phi)
    #
    # [MATLAB] theta_ext_deg = [-fliplr(theta_deg(2:end)), theta_deg, 360 - fliplr(theta_deg(1:end-1))];
    # [MATLAB] phi_ext_deg   = [phi_deg-360, phi_deg, phi_deg+360];
    theta_ext_deg = np.concatenate([
        -theta_deg[1:][::-1],              # top mirror:    -180°, ..., -1°
        theta_deg,                         # main:              0°, ..., 180°
        360.0 - theta_deg[N_theta - 2::-1],  # bottom mirror: 181°, ..., 360°
    ])
    phi_ext_deg = np.concatenate([
        phi_deg - 360.0, phi_deg, phi_deg + 360.0,
    ])

    # Physical theta row index for each row of the extended grid.
    # [MATLAB] theta_phys_idx = [N_theta:-1:2, 1:N_theta, N_theta-1:-1:1];  (1-based)
    theta_phys_idx = np.concatenate([
        np.arange(N_theta - 1, 0, -1),   # top mirror:    N-1 down to 1
        np.arange(N_theta),              # main:              0 to N-1
        np.arange(N_theta - 2, -1, -1), # bottom mirror: N-2 down to 0
    ])  # shape (3N-2,)

    # Phi column offset (in indices) applied to physical phi for mirrored theta rows.
    # Mirror rows shift phi by 180° to reach the back hemisphere.
    # [MATLAB] phi_offset = [repmat(PHI_180_SHIFT, 1, N_theta-1), zeros(1, N_theta), repmat(PHI_180_SHIFT, 1, N_theta-1)];
    phi_offset = np.concatenate([
        np.full(N_theta - 1, PHI_180_SHIFT, dtype=int),  # top mirror
        np.zeros(N_theta, dtype=int),                     # main
        np.full(N_theta - 1, PHI_180_SHIFT, dtype=int),  # bottom mirror
    ])  # shape (3N-2,)

    # |sin(θ)| on the extended grid — solid-angle weight.
    # Absolute value keeps mirror rows (negative θ_ext) positive.
    # [MATLAB] sin_theta_ext = abs(sin(deg2rad(theta_ext_deg)));
    sin_theta_ext = np.abs(np.sin(np.deg2rad(theta_ext_deg)))

    # ── Sparse mask pre-computation ───────────────────────────────────
    # For each directive compute the angular-window mask on the extended grid,
    # then convert the True positions to physical (theta, phi) index pairs.
    # This allows cost_fn to sample power_grid at only the relevant points.
    directive_mask_data = []
    for directive in directives:
        target_theta_deg = directive["theta"]
        target_phi_deg   = directive.get("phi", DEFAULT_TARGET_PHI_DEG)
        sym_width        = directive.get("width", None)
        theta_width_deg  = directive.get("theta_width", sym_width)
        phi_width_deg    = directive.get("phi_width",   sym_width)

        mask = angular_window_mask(
            theta_ext_deg, phi_ext_deg,
            target_theta_deg, target_phi_deg,
            theta_width_deg, phi_width_deg,
        )  # shape (3N-2, 3*N_phi)

        ext_t, ext_p = np.nonzero(mask)   # row / col indices in extended grid

        # Map extended indices → physical grid indices.
        phys_theta = theta_phys_idx[ext_t]
        # ext_p % N_phi gives the base column in the centre phi tile;
        # adding phi_offset then wrapping by N_phi applies the 180° hemisphere shift.
        phys_phi = (ext_p % N_phi + phi_offset[ext_t]) % N_phi
        # [MATLAB] phys_phi = mod(mod(ext_p-1, N_phi) + phi_offset(ext_t), N_phi) + 1;

        sin_w       = sin_theta_ext[ext_t]
        aggregation = directive.get("aggregation", "mean")
        d_type      = directive["type"]
        d_weight    = directive.get("weight", DEFAULT_DIRECTIVE_WEIGHT)

        # Validate here so errors are raised at build time, not inside the hot loop.
        if d_type not in VALID_DIRECTIVE_TYPES:
            raise ValueError(
                f"Unknown directive type '{d_type}'. Valid types: {VALID_DIRECTIVE_TYPES}."
            )
        if aggregation not in VALID_AGGREGATIONS:
            raise ValueError(
                f"Unknown aggregation '{aggregation}'. Valid options: {VALID_AGGREGATIONS}."
            )

        directive_mask_data.append((phys_theta, phys_phi, sin_w,
                                    aggregation, d_type, d_weight))

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
        af_primary = compute_array_factor(weights_complex, element_patterns_stacked)

        # Pre-compute power grid once; reused for every directive this call.
        # For "total" polarisation, power is the incoherent sum of both polarisations.
        if element_patterns_secondary is not None:
            af_secondary = compute_array_factor(weights_complex, element_patterns_secondary)
            power_grid = np.abs(af_primary) ** 2 + np.abs(af_secondary) ** 2
            # [MATLAB] power_grid = abs(af_primary).^2 + abs(af_secondary).^2;
        else:
            power_grid = np.abs(af_primary) ** 2
            # [MATLAB] power_grid = abs(af_primary) .^ 2;

        # Composite cost: weighted sum over all directives.
        # For each directive, sample power_grid at the pre-computed sparse mask indices.
        total_cost = 0.0
        for phys_theta, phys_phi, sin_w, aggregation, d_type, d_weight in directive_mask_data:
            masked_power = power_grid[phys_theta, phys_phi]  # shape (K,), only mask pts
            term = _directive_cost(masked_power, sin_w, d_type, aggregation)
            total_cost += d_weight * term

        return float(total_cost)

    return cost_fn
