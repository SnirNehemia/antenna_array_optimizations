# ══════════════════════════════════════════════════════════════════
# PLOTTER
# All matplotlib visualization for array pattern, weights, and cost history.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend — safe for headless script use
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle, Patch
from pathlib import Path

from src.cost.cost_function import compute_array_factor

# ────────────────────────── CONSTANTS ─────────────────────────────

# Output filenames for each plot type.
POLAR_PLOT_FILENAME      = "pattern_polar.png"
CARTESIAN_PLOT_FILENAME  = "pattern_cartesian.png"
PROJECTION_PLOT_FILENAME = "pattern_2d.png"
WEIGHTS_PLOT_FILENAME    = "weights.png"
COST_HISTORY_FILENAME    = "cost_history.png"
GIF_FILENAME             = "pattern_evolution.gif"

# Maximum number of GIF frames (strides automatically over iterations).
MAX_GIF_FRAMES = 100

# GIF playback speed in frames per second.
GIF_FPS = 10

# Rendering resolution for all saved figures.
PLOT_DPI = 150

# Colors and transparency for directive window overlays.
PEAK_WINDOW_COLOR = "green"
NULL_WINDOW_COLOR = "red"
WINDOW_ALPHA      = 0.15

# Default dynamic range for polar and Cartesian pattern plots.
DEFAULT_DYNAMIC_RANGE_DB = 40


# ────────────────────────── HELPERS ───────────────────────────────

def _nearest_index(grid_deg, target_deg):
    """Return the index of the nearest grid point to a target angle.

    Args:
        grid_deg (np.ndarray): 1D angle grid in degrees.
        target_deg (float): Target angle in degrees.

    Returns:
        int: Index of the nearest grid point.
    """
    return int(np.argmin(np.abs(grid_deg - target_deg)))
    # [MATLAB] [~, idx] = min(abs(grid_deg - target_deg));


def _extract_cut(array_factor_db_grid, theta_deg, phi_deg, output_config):
    """Extract a 1D pattern cut from the 2D dB grid according to output config.

    Supported cut types:
    - ``"theta_cut"``: elevation sweep at fixed phi. Returns ``(theta_deg, db_cut)``.
    - ``"phi_cut"``:   azimuth sweep at fixed theta. Returns ``(phi_deg, db_cut)``.

    Args:
        array_factor_db_grid (np.ndarray): Pre-computed 10·log10(|AF|²),
            shape (N_theta, N_phi). Units: dB.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        output_config (dict): Output settings from config.yaml. Must contain
            ``plot_cut_type`` (str). For ``"theta_cut"``: reads ``plot_phi_deg`` (float).
            For ``"phi_cut"``: reads ``plot_theta_deg`` (float).

    Returns:
        tuple[np.ndarray, np.ndarray, str]: (angle_axis_deg, power_db_cut, xlabel_str)
            - angle_axis_deg: 1D angle grid for the cut axis. Units: degrees.
            - power_db_cut: 1D power values along the cut. Units: dB.
            - xlabel_str: Axis label string for the plot.
    """
    cut_type = output_config["plot_cut_type"]

    if cut_type == "theta_cut":
        fixed_phi_deg = float(output_config.get("plot_phi_deg", 0.0))
        phi_idx       = _nearest_index(phi_deg, fixed_phi_deg)
        power_db_cut  = array_factor_db_grid[:, phi_idx]
        # [MATLAB] power_db_cut = array_factor_db_grid(:, phi_idx);
        return theta_deg, power_db_cut, f"Theta (deg)  [Phi = {fixed_phi_deg:.1f}°]"

    else:  # "phi_cut"
        fixed_theta_deg = float(output_config.get("plot_theta_deg", 90.0))
        theta_idx       = _nearest_index(theta_deg, fixed_theta_deg)
        power_db_cut    = array_factor_db_grid[theta_idx, :]
        # [MATLAB] power_db_cut = array_factor_db_grid(theta_idx, :);
        return phi_deg, power_db_cut, f"Phi (deg)  [Theta = {fixed_theta_deg:.1f}°]"


def _directive_window_bounds(directive, output_config):
    """Compute the angular window [lo, hi] for a directive on the active cut axis.

    For a theta_cut the window is centred on the directive's theta;
    for a phi_cut it is centred on the directive's phi.
    Uses ``theta_width`` / ``phi_width`` when present; falls back to ``width``.

    Args:
        directive (dict): A beam-shaping directive dict with keys ``type``,
            ``theta``, ``phi`` (optional), ``width`` / ``theta_width`` /
            ``phi_width``.
        output_config (dict): Output settings containing ``plot_cut_type``.

    Returns:
        tuple[float, float, float, str]:
            - center_deg: Window centre on the active angle axis. Units: degrees.
            - lo_deg: Lower window bound. Units: degrees.
            - hi_deg: Upper window bound. Units: degrees.
            - color: Overlay color string (green for peak, red for null).
    """
    cut_type  = output_config["plot_cut_type"]
    sym_width = directive.get("width", None)
    theta_hw  = directive.get("theta_width", sym_width) / 2.0
    phi_hw    = directive.get("phi_width",   sym_width) / 2.0
    color     = PEAK_WINDOW_COLOR if directive["type"] == "peak" else NULL_WINDOW_COLOR

    if cut_type == "theta_cut":
        center_deg = directive["theta"]
        hw = theta_hw
    else:
        center_deg = directive.get("phi", 0.0)
        hw = phi_hw

    return center_deg, center_deg - hw, center_deg + hw, color


def _directive_on_cut(directive, output_config):
    """Return (is_visible, is_front_half) for a directive on the active cut.

    A directive is visible only if its target angle is within the directive's
    own half-width of the cut's fixed angle (i.e. the cut plane intersects
    the angular window). Directives that miss the cut are skipped entirely.

    For theta_cut: checks whether the directive's phi is within half_width of
    fixed_phi (front half) or fixed_phi+180° (back half).
    For phi_cut: checks whether the directive's theta is within half_width of
    fixed_theta.

    Args:
        directive (dict): Directive dict with ``theta``, optional ``phi``,
            and ``width`` keys. Units: degrees.
        output_config (dict): Output settings with ``plot_cut_type``,
            ``plot_phi_deg``, ``plot_theta_deg``.

    Returns:
        tuple[bool, bool]: ``(show_on_plot, on_front_half)``.
            ``on_front_half`` is only meaningful when cut_type == ``"theta_cut"``.
    """
    cut_type  = output_config["plot_cut_type"]
    sym_width = directive.get("width", None)
    theta_hw  = directive.get("theta_width", sym_width) / 2.0
    phi_hw    = directive.get("phi_width",   sym_width) / 2.0
    d_phi     = directive.get("phi", 0.0)
    d_theta   = directive["theta"]

    if cut_type == "theta_cut":
        # A theta-cut at fixed phi: visible if the directive's phi window covers the cut.
        fixed_phi = float(output_config.get("plot_phi_deg", 0.0))
        back_phi  = (fixed_phi + 180.0) % 360.0
        delta_front = abs((d_phi - fixed_phi + 180.0) % 360.0 - 180.0)
        delta_back  = abs((d_phi - back_phi  + 180.0) % 360.0 - 180.0)
        on_front = delta_front <= phi_hw
        on_back  = delta_back  <= phi_hw
        # [MATLAB] on_front = abs(mod(d_phi - fixed_phi + 180, 360) - 180) <= phi_hw;
        return (on_front or on_back), on_front
    else:  # phi_cut
        # A phi-cut at fixed theta: visible if the directive's theta window covers the cut.
        fixed_theta = float(output_config.get("plot_theta_deg", 90.0))
        visible     = abs(d_theta - fixed_theta) <= theta_hw
        # [MATLAB] visible = abs(d_theta - fixed_theta) <= theta_hw;
        return visible, True


def _resolve_output_config(output_config, directives):
    """Resolve ``plot_cut_type: "auto"`` to a concrete cut using the first peak directive.

    Returns a shallow copy of ``output_config`` with ``plot_cut_type`` and
    ``plot_phi_deg`` filled in. The original dict is never mutated.

    Args:
        output_config (dict): Output section from config.yaml.
        directives (list[dict]): Beam-shaping directives.

    Returns:
        dict: Resolved output config with ``plot_cut_type`` set to ``"theta_cut"``
            or ``"phi_cut"`` and the corresponding fixed-angle key set.
    """
    cfg = dict(output_config)
    if cfg.get("plot_cut_type") != "auto":
        return cfg
    target = next((d for d in directives if d["type"] == "peak"), None)
    if target is None and directives:
        target = directives[0]
    cfg["plot_cut_type"] = "theta_cut"
    if target is not None:
        cfg["plot_phi_deg"] = float(target.get("phi", 0.0))
    # [MATLAB] % auto: derive fixed_phi from first peak directive
    return cfg


# ────────────────────────── INDIVIDUAL PLOT FUNCTIONS ─────────────

def save_polar_plot(array_factor_db_grid, theta_deg, phi_deg,
                    directives, output_config, dynamic_range_db, output_path):
    """Save a polar radiation pattern plot normalized to 0 dB at the pattern peak.

    For a theta_cut the plot spans -180° to +180°: the front half (phi = plot_phi_deg)
    is placed at positive angles and the back half (phi + 180°) at negative angles,
    giving a full-circle view. For a phi_cut the 0–360° azimuth sweep already fills
    the circle.

    The radial axis is fixed to [-dynamic_range_db, 0] dB.

    Args:
        array_factor_db_grid (np.ndarray): Pre-computed 10·log10(|AF|²),
            shape (N_theta, N_phi). Units: dB.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        directives (list[dict]): Beam-shaping directives for window overlays.
        output_config (dict): Output settings (``plot_cut_type``, ``plot_phi_deg``,
            ``plot_theta_deg``).
        dynamic_range_db (float): Radial axis span. Units: dB.
        output_path (str | Path): Full path for the saved PNG file.
    """
    cut_type = output_config["plot_cut_type"]

    if cut_type == "theta_cut":
        fixed_phi_deg  = float(output_config.get("plot_phi_deg", 0.0))
        fixed_phi_back = (fixed_phi_deg + 180.0) % 360.0
        phi_idx_front  = _nearest_index(phi_deg, fixed_phi_deg)
        phi_idx_back   = _nearest_index(phi_deg, fixed_phi_back)
        front_db       = array_factor_db_grid[:, phi_idx_front]
        back_db        = array_factor_db_grid[:, phi_idx_back]
        # Back half: theta 0→180 at phi+180° reversed to span -180°→0° on the left side.
        angle_full_deg = np.concatenate([-theta_deg[::-1], theta_deg])
        power_full_db  = np.concatenate([back_db[::-1], front_db])
        # [MATLAB] angle_full_deg = [-fliplr(theta_deg), theta_deg];
        # [MATLAB] power_full_db  = [fliplr(back_db), front_db];
        cut_label = f"Phi = {fixed_phi_deg:.1f}°"
    else:  # phi_cut — already full 0–360°
        angle_full_deg, power_full_db, xlabel_str = _extract_cut(
            array_factor_db_grid, theta_deg, phi_deg, output_config
        )
        cut_label     = xlabel_str.split("[")[1].rstrip("]") if "[" in xlabel_str else ""
        fixed_phi_deg = None

    # Normalize so the pattern peak sits at 0 dB on the radial axis.
    power_db_normalized = power_full_db - np.max(power_full_db)
    # [MATLAB] power_db_normalized = power_full_db - max(power_full_db);

    angle_full_rad = np.deg2rad(angle_full_deg)
    # [MATLAB] angle_full_rad = angle_full_deg * pi / 180;

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))

    # Boresight (theta = 0°) at the top; angles increase clockwise.
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    ax.plot(angle_full_rad, power_db_normalized, linewidth=1.5, label="Array pattern")

    # Shade each directive's angular window — only if it intersects this cut.
    for directive in directives:
        visible, on_front = _directive_on_cut(directive, output_config)
        if not visible:
            continue
        center_deg, lo_deg, hi_deg, color = _directive_window_bounds(directive, output_config)
        label_str = f'{directive["type"]} window @ {center_deg:.1f}°'
        if cut_type == "theta_cut":
            if on_front:
                ax.axvspan(np.deg2rad(lo_deg), np.deg2rad(hi_deg),
                           alpha=WINDOW_ALPHA, color=color, label=label_str)
            else:
                ax.axvspan(np.deg2rad(-hi_deg), np.deg2rad(-lo_deg),
                           alpha=WINDOW_ALPHA, color=color, label=label_str)
        else:
            ax.axvspan(np.deg2rad(lo_deg), np.deg2rad(hi_deg),
                       alpha=WINDOW_ALPHA, color=color, label=label_str)
        # [MATLAB] patch(deg2rad([lo lo hi hi]), ylim_vals, color, 'FaceAlpha', WINDOW_ALPHA);

    ax.set_rlim(-dynamic_range_db, 0)
    ax.set_rticks(np.arange(-dynamic_range_db, 1, 10).tolist())
    ax.set_title(f"Array Pattern — {cut_label}", pad=15)
    ax.legend(loc="lower right", fontsize=7)

    fig.savefig(output_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)


def save_cartesian_plot(array_factor_db_grid, theta_deg, phi_deg,
                        directives, output_config, dynamic_range_db, output_path):
    """Save a Cartesian (dB vs angle) pattern plot normalized to 0 dB at peak.

    For a theta_cut the x-axis spans -180° to +180°: the phi+180° back-half cut
    is placed at negative angles. For a phi_cut the 0–360° range is used as-is.

    The y-axis is fixed to [-dynamic_range_db, 0] dB.

    Args:
        array_factor_db_grid (np.ndarray): Pre-computed 10·log10(|AF|²),
            shape (N_theta, N_phi). Units: dB.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        directives (list[dict]): Beam-shaping directives for window overlays.
        output_config (dict): Output settings (``plot_cut_type``, ``plot_phi_deg``,
            ``plot_theta_deg``).
        dynamic_range_db (float): Y-axis span. Units: dB.
        output_path (str | Path): Full path for the saved PNG file.
    """
    cut_type = output_config["plot_cut_type"]

    if cut_type == "theta_cut":
        fixed_phi_deg  = float(output_config.get("plot_phi_deg", 0.0))
        fixed_phi_back = (fixed_phi_deg + 180.0) % 360.0
        phi_idx_front  = _nearest_index(phi_deg, fixed_phi_deg)
        phi_idx_back   = _nearest_index(phi_deg, fixed_phi_back)
        front_db       = array_factor_db_grid[:, phi_idx_front]
        back_db        = array_factor_db_grid[:, phi_idx_back]
        angle_full_deg = np.concatenate([-theta_deg[::-1], theta_deg])
        power_full_db  = np.concatenate([back_db[::-1], front_db])
        xlabel_str = (
            f"Theta (deg)  "
            f"[left: Phi = {fixed_phi_back:.1f}° | right: Phi = {fixed_phi_deg:.1f}°]"
        )
    else:  # phi_cut
        angle_full_deg, power_full_db, xlabel_str = _extract_cut(
            array_factor_db_grid, theta_deg, phi_deg, output_config
        )
        fixed_phi_deg = None

    # Normalize so the pattern peak sits at 0 dB.
    power_db_normalized = power_full_db - np.max(power_full_db)
    # [MATLAB] power_db_normalized = power_full_db - max(power_full_db);

    fig, ax = plt.subplots(figsize=(8, 4))

    ax.plot(angle_full_deg, power_db_normalized, linewidth=1.5, label="Array pattern")

    # Shade directive windows — only if the directive intersects this cut.
    for directive in directives:
        visible, on_front = _directive_on_cut(directive, output_config)
        if not visible:
            continue
        center_deg, lo_deg, hi_deg, color = _directive_window_bounds(directive, output_config)
        label_str = f'{directive["type"]} window @ {center_deg:.1f}°'

        if cut_type == "theta_cut":
            if on_front:
                ax.axvspan(lo_deg, hi_deg, alpha=WINDOW_ALPHA, color=color, label=label_str)
                ax.axvline(center_deg, color=color, linestyle="--", linewidth=0.8, alpha=0.7)
            else:
                ax.axvspan(-hi_deg, -lo_deg, alpha=WINDOW_ALPHA, color=color, label=label_str)
                ax.axvline(-center_deg, color=color, linestyle="--", linewidth=0.8, alpha=0.7)
        else:
            ax.axvspan(lo_deg, hi_deg, alpha=WINDOW_ALPHA, color=color, label=label_str)
            ax.axvline(center_deg, color=color, linestyle="--", linewidth=0.8, alpha=0.7)
        # [MATLAB] xline(center_deg, '--', 'Color', color, 'Alpha', 0.7);

    if cut_type == "theta_cut":
        ax.set_xlim(-180, 180)
    ax.set_ylim(-dynamic_range_db, 0)
    ax.set_xlabel(xlabel_str)
    ax.set_ylabel("Power (dB, normalized to peak)")
    ax.set_title("Array Pattern — Cartesian")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)


def save_2d_projection_plot(array_factor_db_grid, theta_deg, phi_deg,
                             directives, dynamic_range_db, output_path,
                             equal_area=False):
    """Save a 2D theta-phi heatmap of the full radiation pattern with directive markers.

    The full (N_theta × N_phi) grid is normalized to 0 dB at the global peak and
    clamped at -dynamic_range_db. Directives are overlaid as filled rectangles.

    When ``equal_area`` is False (default), the y-axis uses linear θ (0–180°).
    Set ``equal_area=True`` for the cos(θ) equal-area cylindrical projection.

    Args:
        array_factor_db_grid (np.ndarray): Pre-computed 10·log10(|AF|²),
            shape (N_theta, N_phi). Units: dB.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        directives (list[dict]): Beam-shaping directives. Each must have ``type``,
            ``theta``, and optionally ``phi`` (default 0.0), ``theta_width`` /
            ``phi_width`` / ``width``.
        dynamic_range_db (float): Color axis span. Units: dB.
        output_path (str | Path): Full path for the saved PNG file.
        equal_area (bool): Use cos(θ) y-axis for equal solid-angle representation.
            Default: False (linear θ axis).
    """
    # Normalize full 2D grid to peak = 0 dB, then clamp at floor.
    grid_norm = array_factor_db_grid - np.max(array_factor_db_grid)
    grid_norm = np.clip(grid_norm, -dynamic_range_db, 0)
    # [MATLAB] grid_norm = array_factor_db_grid - max(array_factor_db_grid(:));
    # [MATLAB] grid_norm = max(grid_norm, -dynamic_range_db);

    if equal_area:
        # Equal-area cylindrical projection: y = cos(θ).
        # θ=0° (zenith) → cos=+1 (top), θ=90° (horizon) → cos=0, θ=180° → cos=-1 (bottom).
        # [MATLAB] y_coords = cosd(theta_deg);
        y_coords = np.cos(np.deg2rad(theta_deg))
        y_label  = "Elevation θ  [equal-area cos θ]"
        y_lim    = (-1.0, 1.0)
        tick_theta = np.array([0, 30, 60, 90, 120, 150, 180])
        tick_y     = np.cos(np.deg2rad(tick_theta))
        tick_labels = [f"{t}°" for t in tick_theta]
    else:
        y_coords = theta_deg
        y_label  = "Elevation θ (°)"
        y_lim    = (180.0, 0.0)  # inverted: θ=0° (boresight) at top
        tick_theta  = None
        tick_y      = None
        tick_labels = None

    fig, ax = plt.subplots(figsize=(9, 5))

    mesh = ax.pcolormesh(
        phi_deg, y_coords, grid_norm,
        cmap="jet", vmin=-dynamic_range_db, vmax=0, shading="auto",
    )
    # [MATLAB] imagesc(phi_deg, y_coords, grid_norm); colormap jet; caxis([-dynamic_range_db 0]);

    cbar = fig.colorbar(mesh, ax=ax, label="dB (normalized to peak)")
    cbar.set_ticks(np.arange(-dynamic_range_db, 1, 10).tolist())

    if tick_y is not None:
        ax.set_yticks(tick_y)
        ax.set_yticklabels(tick_labels)

    # Overlay each directive's actual optimization window on the heatmap.
    # Physical masks are built with the same extended-grid logic used by the
    # cost function, so wrap-around at φ=0°/360° and pole-crossing at θ=0°/180°
    # are shown exactly as the optimizer sees them.
    from src.cost.cost_function import build_directive_physical_masks  # local import avoids circular dep
    phys_masks = build_directive_physical_masks(theta_deg, phi_deg, directives)

    legend_handles = []
    for directive, phys_mask in zip(directives, phys_masks):
        d_theta   = directive["theta"]
        d_phi     = directive.get("phi", 0.0)
        color     = PEAK_WINDOW_COLOR if directive["type"] == "peak" else NULL_WINDOW_COLOR
        type_str  = "Peak" if directive["type"] == "peak" else "Null"
        label_str = f'{type_str} window  θ={d_theta:.0f}°, φ={d_phi:.0f}°'

        if phys_mask.any():
            mask_f = phys_mask.astype(float)
            ax.contourf(phi_deg, y_coords, mask_f, levels=[0.5, 1.5],
                        colors=[color], alpha=0.30, zorder=4)
            ax.contour(phi_deg, y_coords, mask_f, levels=[0.5],
                       colors=[color], linewidths=2, zorder=5)
            # [MATLAB] ... use imagesc overlay or patch array for mask visualization

        d_y = float(np.cos(np.deg2rad(d_theta))) if equal_area else d_theta
        ax.scatter(d_phi, d_y, marker="+", s=120, color=color, linewidths=2, zorder=6)
        legend_handles.append(
            Patch(facecolor=color, edgecolor=color, alpha=0.30, label=label_str)
        )

    ax.set_xlabel("Phi (deg)")
    ax.set_ylabel(y_label)
    ax.set_title("Array Pattern — 2D Projection")
    ax.legend(handles=legend_handles, fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.2, color="white")

    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(*y_lim)

    fig.tight_layout()
    fig.savefig(output_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)


def save_weight_plots(weights_complex, output_path):
    """Save amplitude and phase bar charts for the optimized complex weights.

    Produces a single figure with two vertically-stacked subplots:
    the top shows element amplitude |w_n|, the bottom shows element phase ∠w_n
    in degrees.

    Args:
        weights_complex (np.ndarray): Optimized complex weights, shape (N_elements,).
            Units: dimensionless (V/V).
        output_path (str | Path): Full path for the saved PNG file.
    """
    n_elements      = len(weights_complex)
    element_indices = np.arange(n_elements)
    # [MATLAB] element_indices = 0 : n_elements - 1;

    amplitudes_linear = np.abs(weights_complex)
    phases_deg        = np.angle(weights_complex, deg=True)
    # [MATLAB] phases_deg = angle(weights_complex) * 180 / pi;

    fig, (ax_amp, ax_phase) = plt.subplots(2, 1, figsize=(8, 5))

    # Amplitude subplot.
    ax_amp.bar(element_indices, amplitudes_linear, width=0.7)
    ax_amp.set_xlabel("Element index")
    ax_amp.set_ylabel("Amplitude |w_n|")
    ax_amp.set_title("Weight Amplitudes")
    ax_amp.set_xticks(element_indices)
    ax_amp.grid(True, axis="y", alpha=0.3)

    # Phase subplot.
    ax_phase.bar(element_indices, phases_deg, width=0.7)
    ax_phase.set_xlabel("Element index")
    ax_phase.set_ylabel("Phase ∠w_n (deg)")
    ax_phase.set_title("Weight Phases")
    ax_phase.set_xticks(element_indices)
    ax_phase.set_ylim(-180, 180)
    ax_phase.axhline(0, color="black", linewidth=0.5)
    ax_phase.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)


def save_cost_history_plot(all_cost_histories, best_run_index, all_run_labels, output_path):
    """Save the optimizer convergence plot showing all restart trajectories.

    Each restart is drawn as a thin grey line. The best run is highlighted as a
    bold blue line labeled with its initialization type.

    Args:
        all_cost_histories (list[list[float]]): Cost history per restart, in order.
        best_run_index (int): Index into all_cost_histories for the best run.
        all_run_labels (list[str]): Human-readable label for each restart.
        output_path (str | Path): Full path for the saved PNG file.
    """
    fig, ax = plt.subplots(figsize=(8, 4))

    first_other = True
    for i, hist in enumerate(all_cost_histories):
        iters = np.arange(len(hist))
        # [MATLAB] iterations = 0 : length(hist) - 1;
        if i == best_run_index:
            ax.plot(iters, hist, linewidth=2, color="C0", zorder=3,
                    label=f"Best — {all_run_labels[i]}")
        else:
            n_other = len(all_cost_histories) - 1
            lbl = f"Other runs ({n_other} total)" if first_other else "_nolegend_"
            ax.plot(iters, hist, linewidth=0.7, color="gray", alpha=0.4,
                    label=lbl, zorder=1)
            first_other = False

    ax.set_xlabel("Iteration")
    ax.set_ylabel("Cost J")
    ax.set_title("Optimizer Convergence — All Restarts")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)


# ────────────────────────── GIF ANIMATION ─────────────────────────

def save_pattern_gif(weights_history, element_patterns_stacked,
                     theta_deg, phi_deg, directives,
                     cost_history, output_config, output_dir):
    """Save an animated GIF showing the radiation pattern evolving over optimizer iterations.

    Each frame renders the 2D heatmap (equal-area projection) with directive overlays
    and a cost-convergence panel. Iteration count is automatically strided to stay
    within ``gif_max_frames``.

    Requires the ``pillow`` package (``pip install pillow``).

    Args:
        weights_history (list[np.ndarray]): Complex weight vectors recorded at each
            optimizer iteration, shape (N_iterations, N_elements). Units: V/V.
        element_patterns_stacked (np.ndarray): Complex element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        directives (list[dict]): Beam-shaping directives for overlays.
        cost_history (list[float]): Cost value per iteration (same length as
            ``weights_history``).
        output_config (dict): Output section from config.yaml. Reads
            ``gif_max_frames`` (int, default 100) and ``plot_dynamic_range_db``.
        output_dir (Path): Directory where the GIF is saved.
    """
    if not weights_history:
        return

    try:
        from PIL import Image  # noqa: F401 — ensure pillow is importable
    except ImportError:
        print("  WARNING: pillow not installed — skipping GIF. Run: pip install pillow")
        return

    max_frames      = int(output_config.get("gif_max_frames", MAX_GIF_FRAMES))
    dynamic_range   = float(output_config.get("plot_dynamic_range_db", DEFAULT_DYNAMIC_RANGE_DB))
    n_iters         = len(weights_history)
    stride          = max(1, n_iters // max_frames)
    frame_indices   = list(range(0, n_iters, stride))
    if (n_iters - 1) not in frame_indices:
        frame_indices.append(n_iters - 1)   # always include the final state

    fig, (ax_map, ax_cost) = plt.subplots(
        1, 2, figsize=(12, 5),
        gridspec_kw={"width_ratios": [2, 1]},
    )

    # ── Initialise heatmap on ax_map ──────────────────────────────────
    dummy_grid = np.zeros((len(theta_deg), len(phi_deg)))
    mesh = ax_map.pcolormesh(
        phi_deg, theta_deg, dummy_grid,
        cmap="jet", vmin=-dynamic_range, vmax=0, shading="auto",
    )
    fig.colorbar(mesh, ax=ax_map, label="dB (normalized to peak)")
    ax_map.set_xlabel("Phi (deg)")
    ax_map.set_ylabel("Elevation θ (°)")
    ax_map.set_xlim(0.0, 360.0)
    ax_map.set_ylim(180.0, 0.0)  # inverted: θ=0° (boresight) at top
    ax_map.grid(True, alpha=0.2, color="white")
    title_text = ax_map.set_title("")

    # Directive overlays — static artists drawn once onto the background.
    # Physical masks match the cost function windows exactly, including
    # phi 0°/360° wrap-around and theta pole-crossing at both poles.
    from src.cost.cost_function import build_directive_physical_masks  # local import
    phys_masks = build_directive_physical_masks(theta_deg, phi_deg, directives)
    for directive, phys_mask in zip(directives, phys_masks):
        d_theta = directive["theta"]
        d_phi   = directive.get("phi", 0.0)
        color   = PEAK_WINDOW_COLOR if directive["type"] == "peak" else NULL_WINDOW_COLOR
        if phys_mask.any():
            mask_f = phys_mask.astype(float)
            ax_map.contourf(phi_deg, theta_deg, mask_f, levels=[0.5, 1.5],
                            colors=[color], alpha=0.30, zorder=4)
            ax_map.contour(phi_deg, theta_deg, mask_f, levels=[0.5],
                           colors=[color], linewidths=2, zorder=5)
        ax_map.scatter(d_phi, d_theta, marker="+", s=100, color=color,
                       linewidths=2, zorder=6)

    # ── Initialise cost plot on ax_cost ───────────────────────────────
    cost_arr = np.array(cost_history)
    iters    = np.arange(len(cost_arr))
    ax_cost.plot(iters, cost_arr, color="gray", alpha=0.5, linewidth=1)
    cost_dot, = ax_cost.plot([], [], "bo", markersize=6)
    ax_cost.set_xlabel("Iteration")
    ax_cost.set_ylabel("Cost J")
    ax_cost.set_title("Convergence")
    ax_cost.grid(True, alpha=0.3)

    fig.tight_layout()

    def _update(frame_idx):
        weights = weights_history[frame_idx]
        # Compute normalised dB grid for this weight vector.
        af         = compute_array_factor(weights, element_patterns_stacked)
        power      = np.abs(af) ** 2
        power_db   = 10.0 * np.log10(np.maximum(power, 1e-30))
        grid_norm  = np.clip(power_db - power_db.max(), -dynamic_range, 0.0)
        mesh.set_array(grid_norm.ravel())
        # [MATLAB] set(mesh_handle, 'CData', grid_norm);

        cost_val = cost_history[frame_idx] if frame_idx < len(cost_history) else float("nan")
        title_text.set_text(f"Iter {frame_idx}   J = {cost_val:.4f}")

        # Cost dot on convergence plot.
        cost_dot.set_data([frame_idx], [cost_val])

        return [mesh, title_text, cost_dot]

    anim = animation.FuncAnimation(
        fig, _update, frames=frame_indices, blit=True, repeat=False,
    )
    interval_ms = int(1000 / GIF_FPS)
    writer = animation.PillowWriter(fps=GIF_FPS)
    anim.save(str(output_dir / GIF_FILENAME), writer=writer)
    plt.close(fig)
    # [MATLAB] MATLAB's imwrite() loop or VideoWriter can produce a similar result.


# ────────────────────────── ORCHESTRATOR ──────────────────────────

def save_all_plots(array_factor_db_grid, theta_deg, phi_deg,
                   weights_complex, directives,
                   all_cost_histories, best_run_index, all_run_labels,
                   output_config, output_dir):
    """Save all enabled plots to the output directory.

    Checks each ``save_*`` flag in ``output_config`` and calls the corresponding
    individual plot function. Files are written to ``output_dir`` under fixed
    filenames defined by the module constants.

    The ``array_factor_db_grid`` must be pre-computed by the caller as:
        ``10 * log10(max(|AF|², 1e-30))``
    to keep the plotter free of array factor computation logic.

    Args:
        array_factor_db_grid (np.ndarray): Pre-computed 10·log10(|AF|²),
            shape (N_theta, N_phi). Units: dB.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        weights_complex (np.ndarray): Optimized complex weights, shape (N_elements,).
            Units: dimensionless (V/V).
        directives (list[dict]): Beam-shaping directives.
        all_cost_histories (list[list[float]]): Cost history per optimizer restart.
        best_run_index (int): Index of the best restart in all_cost_histories.
        all_run_labels (list[str]): Human-readable label for each restart.
        output_config (dict): Output section from config.yaml. Recognized keys:
            ``save_polar_plot``, ``save_cartesian_plot``, ``save_2d_projection_plot``,
            ``save_weight_plots``, ``save_cost_history_plot``, ``plot_cut_type``,
            ``plot_phi_deg``, ``plot_theta_deg``, ``plot_dynamic_range_db``.
        output_dir (Path): Directory where plot files are saved. Must exist.
    """
    output_dir       = Path(output_dir)
    output_config    = _resolve_output_config(output_config, directives)
    dynamic_range_db = float(output_config.get("plot_dynamic_range_db", DEFAULT_DYNAMIC_RANGE_DB))

    if output_config.get("save_polar_plot", False):
        save_polar_plot(
            array_factor_db_grid, theta_deg, phi_deg, directives,
            output_config, dynamic_range_db,
            output_dir / POLAR_PLOT_FILENAME,
        )

    if output_config.get("save_cartesian_plot", False):
        save_cartesian_plot(
            array_factor_db_grid, theta_deg, phi_deg, directives,
            output_config, dynamic_range_db,
            output_dir / CARTESIAN_PLOT_FILENAME,
        )

    if output_config.get("save_2d_projection_plot", False):
        equal_area = bool(output_config.get("plot_equal_area", True))
        save_2d_projection_plot(
            array_factor_db_grid, theta_deg, phi_deg, directives,
            dynamic_range_db,
            output_dir / PROJECTION_PLOT_FILENAME,
            equal_area=equal_area,
        )

    if output_config.get("save_weight_plots", False):
        save_weight_plots(
            weights_complex,
            output_dir / WEIGHTS_PLOT_FILENAME,
        )

    if output_config.get("save_cost_history_plot", False):
        save_cost_history_plot(
            all_cost_histories, best_run_index, all_run_labels,
            output_dir / COST_HISTORY_FILENAME,
        )
