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
from matplotlib.patches import Rectangle
from pathlib import Path

# ────────────────────────── CONSTANTS ─────────────────────────────

# Output filenames for each plot type.
POLAR_PLOT_FILENAME      = "pattern_polar.png"
CARTESIAN_PLOT_FILENAME  = "pattern_cartesian.png"
PROJECTION_PLOT_FILENAME = "pattern_2d.png"
WEIGHTS_PLOT_FILENAME    = "weights.png"
COST_HISTORY_FILENAME    = "cost_history.png"

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

    Args:
        directive (dict): A beam-shaping directive dict with keys ``type``,
            ``theta``, ``phi`` (optional), and ``width``.
        output_config (dict): Output settings containing ``plot_cut_type``.

    Returns:
        tuple[float, float, float, str]:
            - center_deg: Window centre on the active angle axis. Units: degrees.
            - lo_deg: Lower window bound. Units: degrees.
            - hi_deg: Upper window bound. Units: degrees.
            - color: Overlay color string (green for peak, red for null).
    """
    cut_type       = output_config["plot_cut_type"]
    half_width_deg = directive["width"] / 2.0
    color          = PEAK_WINDOW_COLOR if directive["type"] == "peak" else NULL_WINDOW_COLOR

    if cut_type == "theta_cut":
        center_deg = directive["theta"]
    else:
        center_deg = directive.get("phi", 0.0)

    return center_deg, center_deg - half_width_deg, center_deg + half_width_deg, color


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
    cut_type = output_config["plot_cut_type"]
    half_w   = directive["width"] / 2.0
    d_phi    = directive.get("phi", 0.0)
    d_theta  = directive["theta"]

    if cut_type == "theta_cut":
        fixed_phi = float(output_config.get("plot_phi_deg", 0.0))
        back_phi  = (fixed_phi + 180.0) % 360.0
        delta_front = abs((d_phi - fixed_phi + 180.0) % 360.0 - 180.0)
        delta_back  = abs((d_phi - back_phi  + 180.0) % 360.0 - 180.0)
        on_front = delta_front <= half_w
        on_back  = delta_back  <= half_w
        # [MATLAB] on_front = abs(mod(d_phi - fixed_phi + 180, 360) - 180) <= half_w;
        return (on_front or on_back), on_front
    else:  # phi_cut
        fixed_theta = float(output_config.get("plot_theta_deg", 90.0))
        visible     = abs(d_theta - fixed_theta) <= half_w
        # [MATLAB] visible = abs(d_theta - fixed_theta) <= half_w;
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
                             directives, dynamic_range_db, output_path):
    """Save a 2D theta-phi heatmap of the full radiation pattern with directive markers.

    The full (N_theta × N_phi) grid is normalized to 0 dB at the global peak and
    clamped at -dynamic_range_db. Each directive is overlaid as a scatter marker:
    green star (``*``) for peak, red cross (``x``) for null.

    Args:
        array_factor_db_grid (np.ndarray): Pre-computed 10·log10(|AF|²),
            shape (N_theta, N_phi). Units: dB.
        theta_deg (np.ndarray): Elevation angle grid, shape (N_theta,). Units: degrees.
        phi_deg (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: degrees.
        directives (list[dict]): Beam-shaping directives. Each must have ``type``,
            ``theta``, and optionally ``phi`` (default 0.0).
        dynamic_range_db (float): Color axis span. Units: dB.
        output_path (str | Path): Full path for the saved PNG file.
    """
    # Normalize full 2D grid to peak = 0 dB, then clamp at floor.
    grid_norm = array_factor_db_grid - np.max(array_factor_db_grid)
    grid_norm = np.clip(grid_norm, -dynamic_range_db, 0)
    # [MATLAB] grid_norm = array_factor_db_grid - max(array_factor_db_grid(:));
    # [MATLAB] grid_norm = max(grid_norm, -dynamic_range_db);

    fig, ax = plt.subplots(figsize=(9, 5))

    mesh = ax.pcolormesh(
        phi_deg, theta_deg, grid_norm,
        cmap="jet", vmin=-dynamic_range_db, vmax=0, shading="auto",
    )
    # [MATLAB] imagesc(phi_deg, theta_deg, grid_norm); colormap jet; caxis([-dynamic_range_db 0]);

    cbar = fig.colorbar(mesh, ax=ax, label="dB (normalized to peak)")
    cbar.set_ticks(np.arange(-dynamic_range_db, 1, 10).tolist())

    # Overlay each directive as a rectangle showing the full angular window.
    for directive in directives:
        d_theta   = directive["theta"]
        d_phi     = directive.get("phi", 0.0)
        half_w    = directive["width"] / 2.0
        color     = PEAK_WINDOW_COLOR if directive["type"] == "peak" else NULL_WINDOW_COLOR
        type_str  = "Peak" if directive["type"] == "peak" else "Null"
        label_str = f'{type_str} window  θ={d_theta:.0f}°, φ={d_phi:.0f}°'

        rect = Rectangle(
            (d_phi - half_w, d_theta - half_w),
            2 * half_w, 2 * half_w,
            linewidth=2, edgecolor=color, facecolor=color,
            alpha=0.30, zorder=4, label=label_str,
        )
        ax.add_patch(rect)
        # Cross at the exact target center, visible inside the filled rectangle.
        ax.scatter(d_phi, d_theta, marker="+", s=120, color=color,
                   linewidths=2, zorder=5)
        # [MATLAB] rectangle('Position', [d_phi-half_w, d_theta-half_w, 2*half_w, 2*half_w], ...
        # [MATLAB]   'EdgeColor', color, 'FaceColor', color, 'FaceAlpha', 0.30);

    ax.set_xlabel("Phi (deg)")
    ax.set_ylabel("Theta (deg)")
    ax.set_title("Array Pattern — 2D Projection")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.2, color="white")

    # Fix axes to full-sphere coverage regardless of data resolution.
    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(0.0, 180.0)

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
        save_2d_projection_plot(
            array_factor_db_grid, theta_deg, phi_deg, directives,
            dynamic_range_db,
            output_dir / PROJECTION_PLOT_FILENAME,
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
