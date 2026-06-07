# ══════════════════════════════════════════════════════════════════
# OPTIMIZER
# L-BFGS-B gradient descent with optional multi-start for complex weight optimization.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────
import numpy as np
import scipy.optimize

from src.cost.cost_function import build_cost_function, x_to_weights, weights_to_x

# ────────────────────────── CONSTANTS ─────────────────────────────

# Required keys that must be present in the optimizer config dict.
REQUIRED_CONFIG_KEYS = ("max_iterations", "cost_tolerance", "n_restarts")

# Default values for optional optimizer config keys.
DEFAULT_AMPLITUDE_BOUNDS         = None   # no amplitude bounds by default
DEFAULT_PHASE_ONLY               = False
DEFAULT_AMPLITUDE_ONLY           = False
DEFAULT_USE_UNIFORM_INIT         = True   # include uniform-weights starting point
DEFAULT_USE_SINGLE_ELEMENT_INIT  = True   # include per-element single-active starting points
DEFAULT_GRADIENT_TOLERANCE       = 1e-5   # gradient-norm stopping threshold (scipy gtol)

# L-BFGS-B method string for scipy.optimize.minimize.
LBFGSB_METHOD = "L-BFGS-B"


# ────────────────────────── HELPERS ───────────────────────────────

def _build_lbfgsb_bounds(n_elements, amplitude_bounds_config, mode):
    """Build the L-BFGS-B bounds list for the optimization variable vector.

    Bound interpretation depends on the optimization mode:

    - ``"standard"``: each Re(w_n) and Im(w_n) bounded to [-a_max, a_max].
      Note: this approximates the amplitude constraint |w_n| ≤ a_max,
      since the exact constraint is non-convex in the Re/Im parameterization.
    - ``"amplitude_only"``: each amplitude a_n bounded directly to [a_min, a_max].
      This is exact (no approximation), as x = [a_0, ..., a_{N-1}].
    - ``"phase_only"``: no bounds (amplitudes are normalized inside the cost function).

    Args:
        n_elements (int): Number of array elements.
        amplitude_bounds_config (list[float] | None): [a_min, a_max] from config,
            or None for unbounded.
        mode (str): Optimization mode: ``"standard"``, ``"phase_only"``, or
            ``"amplitude_only"``.

    Returns:
        list[tuple[float, float]] | None: Bounds for scipy, or None for unbounded.
    """
    if mode == "phase_only" or amplitude_bounds_config is None:
        return None

    a_min, a_max = amplitude_bounds_config

    if mode == "amplitude_only":
        # x has length N; each variable is an amplitude directly.
        return [(a_min, a_max)] * n_elements
        # [MATLAB] lb = repmat(a_min, n_elements, 1); ub = repmat(a_max, n_elements, 1);

    # Standard mode: x has length 2N; bound each Re and Im component to [-a_max, a_max].
    # a_min is not applied here because Re(w) and Im(w) can be negative.
    return [(-a_max, a_max)] * (2 * n_elements)
    # [MATLAB] lb = repmat(-a_max, 2*n_elements, 1); ub = repmat(a_max, 2*n_elements, 1);


def _uniform_initial_x(n_elements, mode):
    """Build the initial variable vector for uniform weights (amplitude=1, phase=0).

    Args:
        n_elements (int): Number of array elements.
        mode (str): Optimization mode. Determines the length of the returned vector.

    Returns:
        np.ndarray: Initial variable vector. Shape (2N,) for standard/phase-only,
            (N,) for amplitude-only. Units: dimensionless.
    """
    if mode == "amplitude_only":
        # Amplitude-only: initialize all amplitudes to 1.0.
        return np.ones(n_elements, dtype=float)
        # [MATLAB] x0 = ones(n_elements, 1);
    else:
        # Standard / phase-only: uniform complex weights → 2N vector.
        uniform_weights = np.ones(n_elements, dtype=complex)
        return weights_to_x(uniform_weights)
        # [MATLAB] x0 = reshape([ones(n_elements,1), zeros(n_elements,1)]', 1, []);


def _random_initial_x(n_elements, mode, rng):
    """Build a random initial variable vector for multi-start diversification.

    For standard and phase-only modes: unit amplitudes with uniformly random phases.
    For amplitude-only mode: amplitudes drawn uniformly from [0, 1].

    Args:
        n_elements (int): Number of array elements.
        mode (str): Optimization mode.
        rng (np.random.Generator): Random number generator instance.

    Returns:
        np.ndarray: Random initial variable vector. Shape (2N,) for standard/
            phase-only, (N,) for amplitude-only. Units: dimensionless.
    """
    if mode == "amplitude_only":
        # Random amplitudes uniform in [0, 1].
        random_amplitudes = rng.uniform(0.0, 1.0, size=n_elements)
        # [MATLAB] random_amplitudes = rand(n_elements, 1);
        return random_amplitudes
    else:
        # Random unit-amplitude weights with uniformly distributed phases.
        # Unit amplitudes are a neutral starting point that avoids bias toward
        # any particular amplitude distribution.
        random_phases_rad = rng.uniform(0.0, 2.0 * np.pi, size=n_elements)
        # [MATLAB] random_phases_rad = rand(n_elements, 1) * 2 * pi;
        random_weights = np.exp(1j * random_phases_rad)
        return weights_to_x(random_weights)


def _single_element_initial_x(n_elements, element_idx, mode):
    """Build an initial variable vector with only one element active (weight = 1+0j).

    Useful for exploring solutions anchored to a single-element beam. Not
    meaningful for phase_only mode because all amplitudes are normalized to 1
    inside the cost function anyway.

    Args:
        n_elements (int): Number of array elements.
        element_idx (int): Index of the active element (0-based).
        mode (str): Optimization mode.

    Returns:
        np.ndarray: Initial variable vector with only element element_idx active.
    """
    if mode == "amplitude_only":
        x = np.zeros(n_elements, dtype=float)
        x[element_idx] = 1.0
        return x
        # [MATLAB] x0 = zeros(n_elements, 1); x0(element_idx+1) = 1;
    else:  # standard or phase_only
        weights = np.zeros(n_elements, dtype=complex)
        weights[element_idx] = 1.0 + 0j
        return weights_to_x(weights)
        # [MATLAB] x0 = zeros(2*n_elements, 1); x0(2*element_idx+1) = 1;


def _run_single_optimization(cost_fn, x_initial, bounds,
                             max_iterations, cost_tolerance, gradient_tolerance):
    """Run a single L-BFGS-B optimization from a given starting point.

    Tracks the cost value and the raw variable vector at the end of each
    iteration via the callback mechanism.

    Args:
        cost_fn (callable): Scalar cost function ``f(x) -> float``.
        x_initial (np.ndarray): Starting point for the optimization.
        bounds (list[tuple] | None): Variable bounds for L-BFGS-B, or None.
        max_iterations (int): Maximum number of L-BFGS-B iterations.
        cost_tolerance (float): Relative cost improvement threshold (``ftol``):
            stops when |f_{k+1} - f_k| / max(|f_k|, 1) < ftol.
        gradient_tolerance (float): Gradient-norm stopping threshold (``gtol``):
            stops when the infinity norm of the projected gradient < gtol.
            Independent of cost_tolerance; whichever fires first wins.

    Returns:
        tuple[scipy.optimize.OptimizeResult, list[float], list[np.ndarray]]:
            - result: Full scipy optimization result object.
            - cost_history: Cost value recorded at the end of each iteration.
            - xk_history: Raw variable vector ``xk`` at the end of each iteration.
    """
    cost_history = []
    xk_history   = []

    def callback(xk):
        # Record cost and variable state at the current iterate.
        cost_history.append(float(cost_fn(xk)))
        xk_history.append(xk.copy())
    # [MATLAB] options.OutputFcn = @(x, optimValues, state) record_cost(x, optimValues);

    result = scipy.optimize.minimize(
        cost_fn,
        x_initial,
        method=LBFGSB_METHOD,
        bounds=bounds,
        callback=callback,
        options={
            "maxiter": max_iterations,
            "ftol":    cost_tolerance,
            "gtol":    gradient_tolerance,
        },
    )
    # [MATLAB] options = optimoptions('fmincon', 'Algorithm', 'sqp', ...
    # [MATLAB]     'MaxIterations', max_iterations, ...
    # [MATLAB]     'FunctionTolerance', cost_tolerance, ...
    # [MATLAB]     'OptimalityTolerance', gradient_tolerance);
    # [MATLAB] result = fmincon(cost_fn, x_initial, [], [], [], [], lb, ub, [], options);

    return result, cost_history, xk_history


# ────────────────────────── PUBLIC INTERFACE ──────────────────────

def run_optimizer(element_patterns_stacked, theta_deg, phi_deg,
                  directives, optimizer_config,
                  element_patterns_secondary=None):
    """Run multi-start L-BFGS-B optimization to find the best complex element weights.

    Executes ``n_restarts`` optimization runs. The first restart always uses uniform
    initial weights; subsequent restarts use random initialization for diversity.
    Returns the result with the lowest final cost across all restarts.

    Args:
        element_patterns_stacked (np.ndarray): Complex element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        directives (list[dict]): Beam-shaping directives (see
            :func:`~src.cost.cost_function.build_cost_function` for schema).
        optimizer_config (dict): Optimizer settings. Required keys:

            - ``max_iterations`` (int): Maximum L-BFGS-B iterations per restart.
            - ``cost_tolerance`` (float): Convergence tolerance (ftol).
            - ``n_restarts`` (int): Number of optimization runs (first is uniform init).

            Optional keys:

            - ``amplitude_bounds`` (list[float] | None): [a_min, a_max], or None.
            - ``phase_only`` (bool): Fix amplitudes to 1, optimize phases. Default False.
            - ``amplitude_only`` (bool): Fix phases to 0, optimize amplitudes. Default False.
            - ``use_uniform_init`` (bool): Include uniform-weights starting point as run 0.
              Default True. If False, run 0 uses random phases like all others.
            - ``use_single_element_init`` (bool): Add one run per element with only that
              element active. Default True. Disable to reduce total run count.
        element_patterns_secondary (np.ndarray | None): Optional cross-polarisation
            pattern stack, shape (N_elements, N_theta, N_phi). When provided the cost
            function minimises ``|AF_copol|² + |AF_cross|²`` (total field power).
            Forwarded directly to :func:`~src.cost.cost_function.build_cost_function`.
            Default: None.

    Returns:
        dict: Optimization result with keys:

            - ``weights_complex`` (np.ndarray, shape (N_elements,), complex):
              Best optimized complex weights. Units: dimensionless (V/V).
            - ``cost_history`` (list[float]): Cost value per iteration for the best run.
            - ``result`` (scipy.optimize.OptimizeResult): Full scipy result object.

    Raises:
        KeyError: If a required config key is missing.
        ValueError: If both ``phase_only`` and ``amplitude_only`` are True.
    """
    # ── Validate and extract config ────────────────────────────────
    for key in REQUIRED_CONFIG_KEYS:
        if key not in optimizer_config:
            raise KeyError(
                f"Missing required optimizer config key: '{key}'. "
                f"Required keys: {REQUIRED_CONFIG_KEYS}."
            )

    max_iterations           = optimizer_config["max_iterations"]
    cost_tolerance           = optimizer_config["cost_tolerance"]
    n_restarts               = optimizer_config["n_restarts"]
    gradient_tolerance       = optimizer_config.get("gradient_tolerance",       DEFAULT_GRADIENT_TOLERANCE)
    amplitude_bounds         = optimizer_config.get("amplitude_bounds",         DEFAULT_AMPLITUDE_BOUNDS)
    phase_only               = optimizer_config.get("phase_only",               DEFAULT_PHASE_ONLY)
    amplitude_only           = optimizer_config.get("amplitude_only",           DEFAULT_AMPLITUDE_ONLY)
    use_uniform_init         = optimizer_config.get("use_uniform_init",         DEFAULT_USE_UNIFORM_INIT)
    use_single_element_init  = optimizer_config.get("use_single_element_init",  DEFAULT_USE_SINGLE_ELEMENT_INIT)

    if phase_only and amplitude_only:
        raise ValueError(
            "'phase_only' and 'amplitude_only' cannot both be True. "
            "Choose one mode or set both to False for full complex optimization."
        )

    # Resolve the optimization mode string used by cost_function.
    if phase_only:
        mode = "phase_only"
    elif amplitude_only:
        mode = "amplitude_only"
    else:
        mode = "standard"

    n_elements = element_patterns_stacked.shape[0]

    # ── Build cost function and bounds ─────────────────────────────
    cost_fn = build_cost_function(
        element_patterns_stacked, theta_deg, phi_deg, directives, mode=mode,
        element_patterns_secondary=element_patterns_secondary,
    )

    bounds = _build_lbfgsb_bounds(n_elements, amplitude_bounds, mode)

    # ── Multi-start optimization ───────────────────────────────────
    # Use a seeded generator so restarts are reproducible across runs.
    rng = np.random.default_rng(seed=0)
    # [MATLAB] rng(0);  % seed the random number generator

    best_result        = None
    best_cost_history  = []
    best_xk_history    = []
    best_final_cost    = np.inf
    best_run_index     = 0
    all_cost_histories = []
    all_run_labels     = []
    run_index          = 0

    # ── Phase 1: user-configured restarts ─────────────────────────
    for restart_index in range(n_restarts):
        if restart_index == 0 and use_uniform_init:
            # First restart: uniform weights — neutral, deterministic start.
            x_initial = _uniform_initial_x(n_elements, mode)
            run_label = "Uniform init"
        else:
            # Random phases for global search diversity.
            x_initial = _random_initial_x(n_elements, mode, rng)
            run_label = f"Random init {restart_index}"
        # [MATLAB] % restart_index == 0 and use_uniform: uniform; else random

        result, cost_history, xk_history = _run_single_optimization(
            cost_fn, x_initial, bounds, max_iterations, cost_tolerance, gradient_tolerance
        )
        all_cost_histories.append(cost_history)
        all_run_labels.append(run_label)

        if result.fun < best_final_cost:
            best_final_cost   = result.fun
            best_result       = result
            best_cost_history = cost_history
            best_xk_history   = xk_history
            best_run_index    = run_index
        run_index += 1

    # ── Phase 2: single-element initializations ────────────────────
    # One run per element: element k active (weight=1), all others zero.
    if use_single_element_init:
        for element_idx in range(n_elements):
            x_initial = _single_element_initial_x(n_elements, element_idx, mode)
            run_label = f"Element {element_idx} init"
            # [MATLAB] x0 = zeros(2*n_elements, 1); x0(2*element_idx+1) = 1;

            result, cost_history, xk_history = _run_single_optimization(
                cost_fn, x_initial, bounds, max_iterations, cost_tolerance
            )
            all_cost_histories.append(cost_history)
            all_run_labels.append(run_label)

            if result.fun < best_final_cost:
                best_final_cost   = result.fun
                best_result       = result
                best_cost_history = cost_history
                best_xk_history   = xk_history
                best_run_index    = run_index
            run_index += 1

    if best_result is None:
        raise ValueError(
            "No optimization runs were performed. "
            "Set n_restarts >= 1 or enable use_single_element_init."
        )

    # ── Decode best solution ───────────────────────────────────────
    def _decode_x(x):
        """Decode raw variable vector x to complex weights without power-normalisation."""
        if mode == "amplitude_only":
            return x.astype(complex)
        w = x_to_weights(x)
        if mode == "phase_only":
            a = np.abs(w)
            w = w / np.where(a > 0.0, a, 1.0)
        return w

    if mode == "amplitude_only":
        # Amplitude-only: x is a real amplitude vector; Im parts are 0.
        weights_complex = best_result.x.astype(complex)
        # [MATLAB] weights_complex = complex(best_result_x, zeros(size(best_result_x)));
    else:
        weights_complex = x_to_weights(best_result.x)

        if mode == "phase_only":
            # Normalize final weights to unit amplitude for a clean result.
            amplitudes = np.abs(weights_complex)
            safe_amplitudes = np.where(amplitudes > 0.0, amplitudes, 1.0)
            weights_complex = weights_complex / safe_amplitudes
            # [MATLAB] weights_complex = weights_complex ./ max(abs(weights_complex), eps);

    # Decode the per-iteration variable history for the best run into complex weights.
    weights_history = [_decode_x(xk) for xk in best_xk_history]

    # Power-normalise final weights and history so downstream code (metrics,
    # plots, CSV) sees the same scale as the cost function's internal normalisation.
    def _power_normalize(w):
        power = float(np.sum(np.abs(w) ** 2))
        return w / np.sqrt(max(power, 1e-30))

    weights_complex = _power_normalize(weights_complex)
    weights_history = [_power_normalize(w) for w in weights_history]

    # Normalise global phase so element 0 is always real-positive.
    # Global phase has no physical meaning; fixing it makes results reproducible
    # and weight CSVs directly comparable across runs and between Python/MATLAB.
    phase0 = np.angle(weights_complex[0])
    weights_complex = weights_complex * np.exp(-1j * phase0)
    weights_history = [w * np.exp(-1j * np.angle(w[0])) for w in weights_history]
    # [MATLAB] weights_complex = weights_complex * exp(-1j * angle(weights_complex(1)));

    return {
        "weights_complex":    weights_complex,
        "cost_history":       best_cost_history,
        "result":             best_result,
        "all_cost_histories": all_cost_histories,
        "all_run_labels":     all_run_labels,
        "best_run_index":     best_run_index,
        "weights_history":    weights_history,
    }
