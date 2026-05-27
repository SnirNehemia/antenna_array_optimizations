# ══════════════════════════════════════════════════════════════════
# MANUAL_WEIGHTS
# Interactive GUI for manual antenna array weight tuning with live
# 2-D radiation pattern and metric feedback.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")  # [MATLAB] MATLAB has its own figure framework (App Designer)
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import yaml

# Add project root to sys.path so src/ modules are importable
# [MATLAB] MATLAB uses addpath(); no sys.path equivalent
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.io.cst_parser import load_element_patterns
from src.cost.cost_function import compute_array_factor, build_directive_physical_masks
from src.metrics.metrics import evaluate_metrics, _compute_directivity_dbi_grid

# ────────────────────────── CONSTANTS ─────────────────────────────

# dB floor for pattern display (peak sits at 0 dB, floor at -DYNAMIC_RANGE_DB)
DYNAMIC_RANGE_DB = 40

# Guard value for log10 to avoid log(0); same as in cost_function.py
LOG10_EPSILON = 1e-30

# Matplotlib pattern appearance
COLORMAP = "jet"
FIGURE_SIZE = (8, 6)          # inches (width, height)
FIGURE_DPI = 100
GRID_ALPHA = 0.20
GRID_LINEWIDTH = 0.5

# Directive overlay appearance — mirrors plotter.py constants
WINDOW_ALPHA = 0.30
PEAK_OVERLAY_COLOR = "green"
NULL_OVERLAY_COLOR = "red"
SCATTER_MARKER_SIZE = 120
PATCH_LINEWIDTH = 2

# Tkinter colour palette
BG_COLOR = "#f2f2f2"
PANEL_BG = "#e2e2e2"
ENTRY_BG_COLOR = "#ffffff"
ENTRY_ERROR_COLOR = "#ffcccc"

# Font definitions
LABEL_FONT = ("Helvetica", 9)
HEADER_FONT = ("Helvetica", 9, "bold")
METRIC_FONT = ("Courier", 10)

# Entry widget widths (characters) — kept narrow to leave room for inline sliders
ENTRY_WIDTH_AMP = 6
ENTRY_WIDTH_PHASE = 7

# Initial weight state applied at startup
INITIAL_AMPLITUDE_LINEAR = 1.0
INITIAL_PHASE_DEG = 0.0

# Minimum accepted amplitude (non-negative physical constraint)
MIN_AMPLITUDE_LINEAR = 0.0

# CSV column names — must match run_optimization._save_weights_csv output
CSV_COL_ELEMENT_INDEX = "element_index"
CSV_COL_AMPLITUDE = "amplitude"
CSV_COL_PHASE_DEG = "phase_deg"

# Default values for a newly added directive row
DEFAULT_DIRECTIVE_TYPE = "peak"
DEFAULT_DIRECTIVE_THETA_DEG = 0.0
DEFAULT_DIRECTIVE_PHI_DEG = 0.0
DEFAULT_DIRECTIVE_WIDTH_DEG = 5.0
DEFAULT_DIRECTIVE_WEIGHT = 1.0

# Required keys checked before the tuner window opens
REQUIRED_CONFIG_KEYS = ("element_patterns_dir", "directives")

# Scrollable weights panel height (pixels)
WEIGHTS_PANEL_HEIGHT_PX = 240

# Amplitude slider range — linear (V/V)
AMP_SLIDER_MIN = 0.0
AMP_SLIDER_MAX = 2.0

# Phase slider range — degrees
PHASE_SLIDER_MIN_DEG = -180.0
PHASE_SLIDER_MAX_DEG = 180.0

# Right-column panel width — wider to accommodate inline sliders
RIGHT_PANEL_WIDTH_PX = 480

# Supported polarisation modes
# - copol / cross: coherent superposition of the matching complex element pattern
# - total: orthogonal power sum |AF_copol|² + |AF_xpol|² (see _recompute_and_redraw)
POLARIZATION_COPOL = "copol"
POLARIZATION_CROSS = "cross"
POLARIZATION_TOTAL = "total"

# Heatmap display modes
DISPLAY_RELATIVE = "relative"        # 0 dB at live peak, floor at -DYNAMIC_RANGE_DB
DISPLAY_ABSOLUTE = "absolute"        # Absolute directivity in dBi, user-set clim

# Default dBi range used when "absolute" display mode is selected
DEFAULT_DBI_MIN = -40.0
DEFAULT_DBI_MAX =  10.0

# Width (characters) of the dBi min/max Entry widgets in the toolbar
ENTRY_WIDTH_DBI = 6

# Colorbar text shown in each display mode
COLORBAR_LABEL_RELATIVE = "dB (normalized to peak)"
COLORBAR_LABEL_ABSOLUTE = "dBi (absolute directivity)"

# ────────────────────────── CONFIG HELPERS ────────────────────────


def _load_config(config_path: Path) -> dict:
    """Load and return the YAML configuration file.

    Args:
        config_path (Path): Absolute or relative path to config.yaml.

    Returns:
        dict: Parsed configuration.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not config_path.exists():
        raise FileNotFoundError(
            f"Config file not found: {config_path}. "
            "Check the --config argument."
        )
    with open(config_path, "r", encoding="utf-8") as file_handle:
        return yaml.safe_load(file_handle)


def _validate_config(config: dict) -> None:
    """Raise KeyError if any required top-level config key is absent.

    Args:
        config (dict): Parsed configuration dictionary.

    Raises:
        KeyError: With the name of the first missing key.
    """
    for key in REQUIRED_CONFIG_KEYS:
        if key not in config:
            raise KeyError(
                f"Required config key '{key}' is missing from config.yaml. "
                "Check the project schema in config.yaml."
            )


# ────────────────────────── MAIN CLASS ────────────────────────────


class ManualWeightsTuner:
    """Interactive tkinter GUI for manual antenna array weight tuning.

    Presents a live 2-D radiation pattern heatmap that updates whenever the
    user changes element weights or directive targets. All antenna physics
    reuses the existing src/ modules unchanged.

    Args:
        config_path (Path): Path to config.yaml.
    """

    def __init__(self, config_path: Path) -> None:
        """Load data from config, build the UI, and show the initial pattern.

        Args:
            config_path (Path): Path to config.yaml.
        """
        self._config_path = config_path
        self._config = _load_config(config_path)
        _validate_config(self._config)

        # Load all element patterns from directory specified in config
        patterns_dir = Path(self._config["element_patterns_dir"])
        element_patterns = load_element_patterns(patterns_dir)

        # Angular grids are identical across all elements (validated by the parser)
        self._theta_deg: np.ndarray = element_patterns[0]["theta_deg"]  # (N_theta,) deg
        self._phi_deg: np.ndarray = element_patterns[0]["phi_deg"]      # (N_phi,)   deg
        self._n_elements: int = len(element_patterns)

        # Pre-compute stacked element patterns for both polarisation channels.
        # copol: complex co-pol field — E_complex = copol_abs · exp(j · copol_phase)
        # cross: complex cross-pol field — same formula on the cross-pol columns
        # The 'total' display mode combines them in power: |AF_copol|² + |AF_xpol|²
        # (orthogonal channels add in power; see _recompute_and_redraw).
        # [MATLAB] avoid list comprehension; use a for-loop to build a 3-D array
        self._element_patterns_copol: np.ndarray = np.stack(
            [p["E_complex"] for p in element_patterns], axis=0
        )
        self._element_patterns_cross: np.ndarray = np.stack(
            [p["cross_complex"] for p in element_patterns], axis=0
        )

        # Active polarisation flag — switched by the polarisation combobox.
        # Stored as a string rather than a stack reference because 'total' uses
        # both stacks together and cannot be represented by a single one.
        self._active_polarization: str = POLARIZATION_COPOL

        # Ground-truth weight vector: complex (N_elements,), dimensionless V/V
        self._weights_complex: np.ndarray = np.ones(
            self._n_elements, dtype=complex
        )

        # Lists of StringVar objects backing the weight Entry widgets
        self._amp_vars: list[tk.StringVar] = []
        self._phase_vars: list[tk.StringVar] = []

        # Corresponding Entry widgets (kept for bg-colour validation feedback)
        self._amp_entries: list[tk.Entry] = []
        self._phase_entries: list[tk.Entry] = []

        # Slider widgets paired with each entry (amplitude and phase per element)
        self._amp_sliders: list[ttk.Scale] = []
        self._phase_sliders: list[ttk.Scale] = []

        # Guard flag: prevents slider callbacks from recursing during programmatic sync
        self._syncing_weight_display: bool = False

        # Polarisation combobox state (assigned in _build_toolbar)
        self._polarization_var: tk.StringVar

        # Display-mode state (assigned in _build_toolbar):
        # - _display_mode_var:  "relative" or "absolute" (dBi)
        # - _dbi_min_var / _dbi_max_var: colorbar clim entries for absolute mode
        self._display_mode_var: tk.StringVar
        self._dbi_min_var: tk.StringVar
        self._dbi_max_var: tk.StringVar
        self._dbi_min_entry: tk.Entry
        self._dbi_max_entry: tk.Entry

        # Matplotlib artists added for directive overlays (cleared on each redraw)
        self._directive_artists: list = []

        # Directive table rows; each is a plain dict keyed by field name → StringVar
        # plus "_frame" → tk.Frame (used for identity-based removal)
        self._directive_rows: list[dict] = []

        # Absolute dBi directivity grid — always computed in _recompute_and_redraw.
        self._dbi_grid: np.ndarray | None = None

        # Most recently displayed grid (in current display-mode units) — used by hover.
        self._last_display_grid: np.ndarray | None = None

        # Frames set during UI build (assigned in _build_* methods)
        self._directives_inner_frame: tk.Frame
        self._metrics_inner_frame: tk.Frame
        self._ax: object            # matplotlib Axes
        self._mesh: object          # QuadMesh returned by pcolormesh
        self._cbar: object          # Colorbar (label/ticks updated on mode change)
        self._canvas: FigureCanvasTkAgg
        self._hover_text: object    # in-axes text annotation for cursor readout
        self._status_label: tk.Label
        self._label_cost: tk.Label
        self._label_peak: tk.Label
        self._label_peak_angle: tk.Label
        self._label_hpbw: tk.Label

        # Build the window
        self._root = tk.Tk()
        self._root.title(f"Manual Pattern Tuner — {config_path.name}")
        self._root.configure(bg=BG_COLOR)
        self._root.geometry("1400x860")
        self._root.minsize(1000, 640)

        self._build_ui()

        # Pre-populate directives from config
        for directive_dict in self._config.get("directives", []):
            self._add_directive_row(directive_dict)

        # Render initial pattern with uniform weights
        self._recompute_and_redraw()

        self._set_status(
            f"{self._n_elements} elements loaded · "
            f"{len(self._directive_rows)} directives active"
        )

    # ── UI CONSTRUCTION ───────────────────────────────────────────

    def _build_ui(self) -> None:
        """Assemble the top-level window: toolbar, main columns, status bar."""
        self._build_toolbar()

        main_frame = tk.Frame(self._root, bg=BG_COLOR)
        main_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=4, pady=4)

        # Left column: matplotlib pattern (~60 % of width, expands)
        left_frame = tk.Frame(main_frame, bg=BG_COLOR)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Right column: controls panel (fixed width)
        right_frame = tk.Frame(main_frame, bg=BG_COLOR, width=RIGHT_PANEL_WIDTH_PX)
        right_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=(4, 0))
        right_frame.pack_propagate(False)

        self._build_pattern_panel(left_frame)
        self._build_weights_panel(right_frame)
        self._build_directives_panel(right_frame)
        self._build_metrics_panel(right_frame)

    def _build_toolbar(self) -> None:
        """Build the top toolbar: config label, polarisation selector, action buttons."""
        toolbar = tk.Frame(self._root, bg=PANEL_BG, relief=tk.FLAT, bd=1)
        toolbar.pack(side=tk.TOP, fill=tk.X)

        tk.Label(
            toolbar,
            text=f"  Config: {self._config_path}",
            font=LABEL_FONT,
            bg=PANEL_BG,
            fg="#555555",
        ).pack(side=tk.LEFT, padx=4, pady=5)

        # Load CSV and Uniform buttons on the right
        ttk.Button(
            toolbar, text="Load Weights CSV", command=self._on_load_csv
        ).pack(side=tk.RIGHT, padx=6, pady=4)

        ttk.Button(
            toolbar, text="Uniform Weights", command=self._on_uniform
        ).pack(side=tk.RIGHT, padx=2, pady=4)

        ttk.Separator(toolbar, orient="vertical").pack(
            side=tk.RIGHT, fill=tk.Y, padx=10, pady=4
        )

        # Polarisation selector — controls which element pattern stack is active
        self._polarization_var = tk.StringVar(value=POLARIZATION_COPOL)
        pol_combo = ttk.Combobox(
            toolbar,
            textvariable=self._polarization_var,
            values=[POLARIZATION_COPOL, POLARIZATION_CROSS, POLARIZATION_TOTAL],
            width=6,
            state="readonly",
        )
        pol_combo.pack(side=tk.RIGHT, pady=4)
        pol_combo.bind(
            "<<ComboboxSelected>>", lambda _: self._on_polarization_change()
        )
        tk.Label(
            toolbar, text="Polarisation:", font=LABEL_FONT, bg=PANEL_BG
        ).pack(side=tk.RIGHT, padx=(0, 2), pady=5)

        ttk.Separator(toolbar, orient="vertical").pack(
            side=tk.RIGHT, fill=tk.Y, padx=10, pady=4
        )

        # ── Display-mode controls ────────────────────────────────────
        # max (dBi) entry — packed first because tk.RIGHT stacks right-to-left.
        self._dbi_max_var = tk.StringVar(value=f"{DEFAULT_DBI_MAX:.1f}")
        self._dbi_max_entry = tk.Entry(
            toolbar, textvariable=self._dbi_max_var,
            width=ENTRY_WIDTH_DBI, font=LABEL_FONT, bg=ENTRY_BG_COLOR,
        )
        self._dbi_max_entry.pack(side=tk.RIGHT, pady=4)
        self._dbi_max_entry.bind("<Return>",   lambda _: self._on_display_mode_change())
        self._dbi_max_entry.bind("<FocusOut>", lambda _: self._on_display_mode_change())
        tk.Label(
            toolbar, text="max (dBi):", font=LABEL_FONT, bg=PANEL_BG
        ).pack(side=tk.RIGHT, padx=(4, 2), pady=5)

        self._dbi_min_var = tk.StringVar(value=f"{DEFAULT_DBI_MIN:.1f}")
        self._dbi_min_entry = tk.Entry(
            toolbar, textvariable=self._dbi_min_var,
            width=ENTRY_WIDTH_DBI, font=LABEL_FONT, bg=ENTRY_BG_COLOR,
        )
        self._dbi_min_entry.pack(side=tk.RIGHT, pady=4)
        self._dbi_min_entry.bind("<Return>",   lambda _: self._on_display_mode_change())
        self._dbi_min_entry.bind("<FocusOut>", lambda _: self._on_display_mode_change())
        tk.Label(
            toolbar, text="min (dBi):", font=LABEL_FONT, bg=PANEL_BG
        ).pack(side=tk.RIGHT, padx=(4, 2), pady=5)

        # Display-mode combobox: relative-to-peak vs absolute dBi
        self._display_mode_var = tk.StringVar(value=DISPLAY_RELATIVE)
        disp_combo = ttk.Combobox(
            toolbar,
            textvariable=self._display_mode_var,
            values=[DISPLAY_RELATIVE, DISPLAY_ABSOLUTE],
            width=8,
            state="readonly",
        )
        disp_combo.pack(side=tk.RIGHT, pady=4)
        disp_combo.bind(
            "<<ComboboxSelected>>", lambda _: self._on_display_mode_change()
        )
        tk.Label(
            toolbar, text="Display:", font=LABEL_FONT, bg=PANEL_BG
        ).pack(side=tk.RIGHT, padx=(0, 2), pady=5)

    def _build_pattern_panel(self, parent: tk.Frame) -> None:
        """Embed the matplotlib 2-D heatmap figure in the left column.

        Args:
            parent (tk.Frame): Left-column frame to pack into.
        """
        lf = ttk.LabelFrame(parent, text="2-D Radiation Pattern")
        lf.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        fig = Figure(figsize=FIGURE_SIZE, dpi=FIGURE_DPI)
        fig.patch.set_facecolor(BG_COLOR)

        self._ax = fig.add_subplot(111)
        self._ax.set_facecolor("black")

        # Initialise mesh with a silent (uniform low) pattern so the axes are
        # fully configured before the first real computation
        n_theta = len(self._theta_deg)
        n_phi = len(self._phi_deg)
        init_grid = np.full((n_theta, n_phi), -DYNAMIC_RANGE_DB, dtype=float)

        # 2-D heatmap: azimuth φ on x-axis, elevation θ on y-axis (0°–180°).
        # [MATLAB] use imagesc or pcolor; pcolormesh has no direct equivalent
        self._mesh = self._ax.pcolormesh(
            self._phi_deg,
            self._theta_deg,
            init_grid,
            cmap=COLORMAP,
            vmin=-DYNAMIC_RANGE_DB,
            vmax=0.0,
            shading="auto",
        )

        # Keep a handle on the colorbar — its label and ticks are updated when
        # the user switches between relative-to-peak and absolute-dBi modes.
        self._cbar = fig.colorbar(
            self._mesh, ax=self._ax, label=COLORBAR_LABEL_RELATIVE
        )
        cbar_ticks = list(range(-DYNAMIC_RANGE_DB, 1, 10))
        self._cbar.set_ticks(cbar_ticks)

        self._ax.set_xlabel("Azimuth  φ (°)")
        self._ax.set_ylabel("Elevation θ (°)")
        self._ax.set_title("Array Factor — manual weights", fontsize=10)
        self._ax.grid(
            True, alpha=GRID_ALPHA, color="white", linewidth=GRID_LINEWIDTH
        )

        # Fix axes to full-sphere coverage regardless of data resolution.
        # θ is inverted so 0° (boresight / zenith) is at the top.
        self._ax.set_xlim(0.0, 360.0)
        self._ax.set_ylim(180.0, 0.0)

        # Cursor readout annotation — overlaid in the top-left of the axes.
        # Invisible until the mouse enters the axes for the first time.
        self._hover_text = self._ax.text(
            0.01, 0.98, "",
            transform=self._ax.transAxes,
            va="top", ha="left",
            fontsize=9, color="white",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.55),
            zorder=10,
            visible=False,
        )

        fig.tight_layout()

        self._canvas = FigureCanvasTkAgg(fig, master=lf)
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Status bar sits directly below the canvas, inside the pattern LabelFrame.
        self._status_label = tk.Label(
            lf, text="", font=LABEL_FONT, bg=PANEL_BG, fg="#555555",
            anchor="w", relief=tk.SUNKEN, bd=1,
        )
        self._status_label.pack(side=tk.BOTTOM, fill=tk.X)

        # Connect mouse events: hover readout + leave to hide the annotation.
        self._canvas.mpl_connect("motion_notify_event", self._on_mouse_move)
        self._canvas.mpl_connect("axes_leave_event",    self._on_mouse_leave)

    def _build_weights_panel(self, parent: tk.Frame) -> None:
        """Build the scrollable element-weights panel.

        Each row shows: element index, amplitude entry, phase entry, Solo button.

        Args:
            parent (tk.Frame): Right-column frame to pack into.
        """
        lf = ttk.LabelFrame(parent, text="Element Weights")
        lf.pack(fill=tk.X, padx=4, pady=(4, 2))

        # Column header
        header = tk.Frame(lf, bg=PANEL_BG)
        header.pack(fill=tk.X, padx=4, pady=(2, 0))
        for col_text, col_width in [
            ("Elem", 5), ("Amplitude", 10), ("Phase (°)", 10)
        ]:
            tk.Label(
                header,
                text=col_text,
                font=HEADER_FONT,
                bg=PANEL_BG,
                width=col_width,
                anchor="center",
            ).pack(side=tk.LEFT)

        ttk.Separator(lf, orient="horizontal").pack(fill=tk.X, padx=4, pady=2)

        # Scrollable content area
        scroll_outer = tk.Frame(lf)
        scroll_outer.pack(fill=tk.BOTH)

        scroll_canvas = tk.Canvas(
            scroll_outer,
            bg=BG_COLOR,
            highlightthickness=0,
            height=WEIGHTS_PANEL_HEIGHT_PX,
        )
        scrollbar = ttk.Scrollbar(
            scroll_outer, orient="vertical", command=scroll_canvas.yview
        )
        scroll_frame = tk.Frame(scroll_canvas, bg=BG_COLOR)

        # Resize scroll region whenever the inner frame changes size
        # [MATLAB] MATLAB App Designer uses fixed-size panels; no dynamic scroll region
        scroll_frame.bind(
            "<Configure>",
            lambda e: scroll_canvas.configure(
                scrollregion=scroll_canvas.bbox("all")
            ),
        )
        scroll_canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        scroll_canvas.configure(yscrollcommand=scrollbar.set)

        scroll_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Windows mouse-wheel scrolling (delta is a multiple of 120)
        # [MATLAB] no direct equivalent; MATLAB figures do not scroll this way
        scroll_canvas.bind(
            "<MouseWheel>",
            lambda e: scroll_canvas.yview_scroll(
                int(-1 * (e.delta / 120)), "units"
            ),
        )

        # One row per element
        for n in range(self._n_elements):
            self._build_weight_row(scroll_frame, n)

    def _build_weight_row(self, parent: tk.Frame, elem_idx: int) -> None:
        """Build one element-weight row: index label, entry+slider pairs, Solo button.

        Each of amplitude and phase gets a narrow text entry for precise numeric
        input paired with a horizontal slider for fast sweeping. The two controls
        are kept in sync bidirectionally.

        Args:
            parent (tk.Frame): Scrollable inner frame.
            elem_idx (int): Zero-based element index. Units: dimensionless.
        """
        row = tk.Frame(parent, bg=BG_COLOR)
        row.pack(fill=tk.X, padx=4, pady=1)

        # Element index label
        tk.Label(
            row, text=f"{elem_idx:>3}", font=LABEL_FONT, bg=BG_COLOR,
            width=4, anchor="e",
        ).pack(side=tk.LEFT, padx=(0, 2))

        # ── Amplitude sub-frame: entry on left, slider fills remaining space ──
        amp_frame = tk.Frame(row, bg=BG_COLOR)
        amp_frame.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 4))

        amp_var = tk.StringVar(value=f"{INITIAL_AMPLITUDE_LINEAR:.4f}")
        amp_entry = tk.Entry(
            amp_frame, textvariable=amp_var, width=ENTRY_WIDTH_AMP,
            font=LABEL_FONT, bg=ENTRY_BG_COLOR,
        )
        amp_entry.pack(side=tk.LEFT)

        # [MATLAB] use uicontrol('Style','slider') for a slider widget
        amp_slider = ttk.Scale(
            amp_frame, from_=AMP_SLIDER_MIN, to=AMP_SLIDER_MAX,
            orient="horizontal",
            command=lambda val, idx=elem_idx: self._on_amp_slider_change(idx, val),
        )
        amp_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(2, 0))

        # ── Phase sub-frame: entry on left, slider fills remaining space ──
        phase_frame = tk.Frame(row, bg=BG_COLOR)
        phase_frame.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 4))

        phase_var = tk.StringVar(value=f"{INITIAL_PHASE_DEG:.2f}")
        phase_entry = tk.Entry(
            phase_frame, textvariable=phase_var, width=ENTRY_WIDTH_PHASE,
            font=LABEL_FONT, bg=ENTRY_BG_COLOR,
        )
        phase_entry.pack(side=tk.LEFT)

        phase_slider = ttk.Scale(
            phase_frame, from_=PHASE_SLIDER_MIN_DEG, to=PHASE_SLIDER_MAX_DEG,
            orient="horizontal",
            command=lambda val, idx=elem_idx: self._on_phase_slider_change(idx, val),
        )
        phase_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(2, 0))

        # Commit entry value on Return or focus loss
        # [MATLAB] use uicontrol Callback property; no FocusOut equivalent
        commit_fn = lambda event, idx=elem_idx: self._on_weight_entry_commit(idx)
        amp_entry.bind("<Return>", commit_fn)
        amp_entry.bind("<FocusOut>", commit_fn)
        phase_entry.bind("<Return>", commit_fn)
        phase_entry.bind("<FocusOut>", commit_fn)

        # Solo button: activate only this element, zero all others
        ttk.Button(
            row, text="Solo",
            command=lambda idx=elem_idx: self._on_solo_element(idx),
            width=4,
        ).pack(side=tk.LEFT)

        # Register all widgets before calling set() — ttk.Scale.set() fires the
        # command callback synchronously on Windows, so the lists must already
        # contain this element's entries when the callback runs.
        self._amp_vars.append(amp_var)
        self._phase_vars.append(phase_var)
        self._amp_entries.append(amp_entry)
        self._phase_entries.append(phase_entry)
        self._amp_sliders.append(amp_slider)
        self._phase_sliders.append(phase_slider)

        # Guard prevents the init set() calls from triggering a recompute before
        # the rest of the UI (pattern panel, metrics labels) is fully built.
        self._syncing_weight_display = True
        amp_slider.set(INITIAL_AMPLITUDE_LINEAR)
        phase_slider.set(INITIAL_PHASE_DEG)
        self._syncing_weight_display = False

    def _build_directives_panel(self, parent: tk.Frame) -> None:
        """Build the directives table with add/remove capability.

        Args:
            parent (tk.Frame): Right-column frame to pack into.
        """
        lf = ttk.LabelFrame(parent, text="Directives")
        lf.pack(fill=tk.X, padx=4, pady=2)

        # Table column header
        header = tk.Frame(lf, bg=PANEL_BG)
        header.pack(fill=tk.X, padx=4, pady=(2, 0))
        col_specs = [
            ("Type", 5), ("θ (°)", 5), ("φ (°)", 5), ("θW(°)", 5), ("φW(°)", 5),
            ("wt", 4), ("", 2), ("Result", 14),
        ]
        for col_text, col_width in col_specs:
            tk.Label(
                header,
                text=col_text,
                font=HEADER_FONT,
                bg=PANEL_BG,
                width=col_width,
                anchor="center",
            ).pack(side=tk.LEFT, padx=1)

        ttk.Separator(lf, orient="horizontal").pack(fill=tk.X, padx=4, pady=2)

        # Container frame for dynamically added rows
        self._directives_inner_frame = tk.Frame(lf, bg=BG_COLOR)
        self._directives_inner_frame.pack(fill=tk.X, padx=4)

        ttk.Button(
            lf,
            text="＋  Add Directive",
            command=lambda: self._add_directive_row(None),
        ).pack(pady=4)

    def _build_metrics_panel(self, parent: tk.Frame) -> None:
        """Build the live metrics display panel.

        Args:
            parent (tk.Frame): Right-column frame to pack into.
        """
        lf = ttk.LabelFrame(parent, text="Metrics")
        lf.pack(fill=tk.X, padx=4, pady=(2, 4))

        self._metrics_inner_frame = tk.Frame(lf, bg=BG_COLOR)
        self._metrics_inner_frame.pack(fill=tk.X, padx=6, pady=4)

        # Fixed metric rows (always visible)
        self._label_cost        = self._make_metric_label("Total J:",     "—")
        self._label_peak        = self._make_metric_label("Global peak:", "—")
        self._label_peak_angle  = self._make_metric_label("Peak angle:",  "—")
        self._label_hpbw        = self._make_metric_label("3 dB HPBW:",  "—")

    def _make_metric_label(self, title: str, initial_value: str) -> tk.Label:
        """Create a key = value row inside the metrics frame.

        Args:
            title (str): Left-side descriptive label.
            initial_value (str): Starting right-side text.

        Returns:
            tk.Label: The value label (updated during live operation).
        """
        row = tk.Frame(self._metrics_inner_frame, bg=BG_COLOR)
        row.pack(fill=tk.X, pady=1)

        tk.Label(
            row,
            text=title,
            font=LABEL_FONT,
            bg=BG_COLOR,
            width=14,
            anchor="w",
        ).pack(side=tk.LEFT)

        value_label = tk.Label(
            row,
            text=initial_value,
            font=METRIC_FONT,
            bg=BG_COLOR,
            anchor="w",
            fg="#222222",
        )
        value_label.pack(side=tk.LEFT)
        return value_label

    # ── DIRECTIVE TABLE MANAGEMENT ────────────────────────────────

    def _add_directive_row(self, directive_dict: dict | None) -> None:
        """Append a new row to the directives table and trigger a redraw.

        Args:
            directive_dict (dict | None): Pre-populated values from config, or
                None to use project defaults. Expected keys: type, theta (deg),
                phi (deg), width (deg), weight.
        """
        # Fall back to default values for any missing field.
        # theta_width / phi_width take precedence over the symmetric "width" shorthand.
        defaults = {
            "type":        DEFAULT_DIRECTIVE_TYPE,
            "theta":       DEFAULT_DIRECTIVE_THETA_DEG,
            "phi":         DEFAULT_DIRECTIVE_PHI_DEG,
            "theta_width": DEFAULT_DIRECTIVE_WIDTH_DEG,
            "phi_width":   DEFAULT_DIRECTIVE_WIDTH_DEG,
            "weight":      DEFAULT_DIRECTIVE_WEIGHT,
        }
        d = directive_dict if directive_dict is not None else {}
        sym_w = d.get("width", None)   # symmetric fallback if per-axis widths absent
        theta_width_default = d.get("theta_width", sym_w if sym_w is not None else defaults["theta_width"])
        phi_width_default   = d.get("phi_width",   sym_w if sym_w is not None else defaults["phi_width"])

        row_frame = tk.Frame(self._directives_inner_frame, bg=BG_COLOR)
        row_frame.pack(fill=tk.X, pady=1)

        # Type combobox (readonly; only "peak" and "null" are valid types)
        type_var = tk.StringVar(value=str(d.get("type", defaults["type"])))
        type_combo = ttk.Combobox(
            row_frame,
            textvariable=type_var,
            values=["peak", "null"],
            width=5,
            state="readonly",
        )
        type_combo.pack(side=tk.LEFT, padx=1)
        type_combo.bind(
            "<<ComboboxSelected>>", lambda _: self._on_directive_change()
        )

        # Numeric fields: theta, phi, theta_width, phi_width, weight
        numeric_field_specs = [
            ("theta",       d.get("theta",   defaults["theta"]),   5),
            ("phi",         d.get("phi",     defaults["phi"]),     5),
            ("theta_width", theta_width_default,                   5),
            ("phi_width",   phi_width_default,                     5),
            ("weight",      d.get("weight",  defaults["weight"]),  4),
        ]
        vars_dict: dict = {"type_var": type_var, "_frame": row_frame}

        for field_name, default_val, entry_width in numeric_field_specs:
            field_var = tk.StringVar(value=f"{float(default_val):.2f}")
            field_entry = tk.Entry(
                row_frame,
                textvariable=field_var,
                width=entry_width,
                font=LABEL_FONT,
                bg=ENTRY_BG_COLOR,
            )
            field_entry.pack(side=tk.LEFT, padx=1)
            field_entry.bind("<Return>",   lambda e: self._on_directive_change())
            field_entry.bind("<FocusOut>", lambda e: self._on_directive_change())
            vars_dict[field_name] = field_var

        # Remove button: destroys this row and removes it from the list
        ttk.Button(
            row_frame,
            text="×",
            command=lambda rf=row_frame: self._remove_directive_row(rf),
            width=2,
        ).pack(side=tk.LEFT, padx=(4, 0))

        # Inline live metric readout (gain in dBi, or gain + null depth for nulls)
        metric_label = tk.Label(
            row_frame, text="—", font=METRIC_FONT, bg=BG_COLOR,
            anchor="w", width=14,
        )
        metric_label.pack(side=tk.LEFT, padx=(4, 0))
        vars_dict["metric_label"] = metric_label

        self._directive_rows.append(vars_dict)
        self._on_directive_change()
        self._update_status_count()

    def _remove_directive_row(self, row_frame: tk.Frame) -> None:
        """Remove a directive row identified by its frame widget.

        Args:
            row_frame (tk.Frame): The container frame of the row to remove.
        """
        # Filter by object identity rather than index (indices shift after removal)
        self._directive_rows = [
            r for r in self._directive_rows if r["_frame"] is not row_frame
        ]
        row_frame.destroy()
        self._on_directive_change()
        self._update_status_count()

    # ── EVENT HANDLERS ────────────────────────────────────────────

    def _on_weight_entry_commit(self, elem_idx: int) -> None:
        """Validate and apply a weight entry edit for element elem_idx.

        Args:
            elem_idx (int): Zero-based index of the edited element.
        """
        amp_text   = self._amp_vars[elem_idx].get().strip()
        phase_text = self._phase_vars[elem_idx].get().strip()

        try:
            amplitude_linear = float(amp_text)
            if amplitude_linear < MIN_AMPLITUDE_LINEAR:
                raise ValueError("Amplitude must be non-negative.")
            phase_deg = float(phase_text)
        except ValueError:
            self._amp_entries[elem_idx].config(bg=ENTRY_ERROR_COLOR)
            self._phase_entries[elem_idx].config(bg=ENTRY_ERROR_COLOR)
            self._set_status(
                f"Element {elem_idx}: invalid entry — "
                "amplitude must be ≥ 0, phase must be a number."
            )
            return

        self._amp_entries[elem_idx].config(bg=ENTRY_BG_COLOR)
        self._phase_entries[elem_idx].config(bg=ENTRY_BG_COLOR)

        # Encode amplitude and phase (degrees) as a complex weight V/V.
        # Internal convention: phase always in radians for computation.
        # [MATLAB] use complex(amplitude * cos(phase_rad), amplitude * sin(phase_rad))
        phase_rad = phase_deg * (np.pi / 180.0)
        self._weights_complex[elem_idx] = amplitude_linear * np.exp(1j * phase_rad)

        # Sync slider positions to the committed entry values.
        # Guard prevents the slider callbacks from firing a second recompute.
        self._syncing_weight_display = True
        self._amp_sliders[elem_idx].set(amplitude_linear)
        self._phase_sliders[elem_idx].set(phase_deg)
        self._syncing_weight_display = False

        self._recompute_and_redraw()

    def _on_solo_element(self, elem_idx: int) -> None:
        """Set element elem_idx to unit weight and zero all other elements.

        Args:
            elem_idx (int): Index of the element to activate alone.
        """
        self._weights_complex = np.zeros(self._n_elements, dtype=complex)
        self._weights_complex[elem_idx] = 1.0 + 0.0j
        self._sync_entries_to_weights()
        self._recompute_and_redraw()

    def _on_uniform(self) -> None:
        """Set all elements to amplitude = 1.0 and phase = 0.0°."""
        self._weights_complex = np.ones(self._n_elements, dtype=complex)
        self._sync_entries_to_weights()
        self._recompute_and_redraw()

    def _on_load_csv(self) -> None:
        """Open a file dialog, parse a weights CSV, and apply the result."""
        filepath = filedialog.askopenfilename(
            title="Load Weights CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not filepath:
            return  # User cancelled

        try:
            loaded_weights = self._parse_weights_csv(filepath)
        except (ValueError, KeyError, OSError) as exc:
            messagebox.showerror("Load Error", str(exc))
            return

        self._weights_complex = loaded_weights
        self._sync_entries_to_weights()
        self._recompute_and_redraw()
        self._set_status(f"Weights loaded from: {filepath}")

    def _on_directive_change(self) -> None:
        """Recompute and redraw after any change to the directives table."""
        self._recompute_and_redraw()

    def _on_amp_slider_change(self, elem_idx: int, value: str) -> None:
        """Apply an amplitude slider drag for element elem_idx.

        Skipped during programmatic sync to prevent recursive callbacks.

        Args:
            elem_idx (int): Zero-based element index.
            value (str): New slider value as a string (provided by ttk.Scale command).
                Units: dimensionless linear amplitude (V/V).
        """
        if self._syncing_weight_display:
            return

        amplitude_linear = max(MIN_AMPLITUDE_LINEAR, round(float(value), 4))
        self._amp_vars[elem_idx].set(f"{amplitude_linear:.4f}")
        self._amp_entries[elem_idx].config(bg=ENTRY_BG_COLOR)

        # Re-encode as complex weight using the current phase from the entry.
        # Parse the entry directly; it was last set by either slider or commit.
        try:
            phase_deg = float(self._phase_vars[elem_idx].get())
        except ValueError:
            phase_deg = 0.0
        phase_rad = phase_deg * (np.pi / 180.0)
        self._weights_complex[elem_idx] = amplitude_linear * np.exp(1j * phase_rad)

        self._recompute_and_redraw()

    def _on_phase_slider_change(self, elem_idx: int, value: str) -> None:
        """Apply a phase slider drag for element elem_idx.

        Skipped during programmatic sync to prevent recursive callbacks.

        Args:
            elem_idx (int): Zero-based element index.
            value (str): New slider value as a string (provided by ttk.Scale command).
                Units: degrees.
        """
        if self._syncing_weight_display:
            return

        phase_deg = round(float(value), 1)
        self._phase_vars[elem_idx].set(f"{phase_deg:.2f}")
        self._phase_entries[elem_idx].config(bg=ENTRY_BG_COLOR)

        # Re-encode as complex weight using the current amplitude from the entry.
        try:
            amplitude_linear = max(
                MIN_AMPLITUDE_LINEAR, float(self._amp_vars[elem_idx].get())
            )
        except ValueError:
            amplitude_linear = 0.0
        phase_rad = phase_deg * (np.pi / 180.0)
        self._weights_complex[elem_idx] = amplitude_linear * np.exp(1j * phase_rad)

        self._recompute_and_redraw()

    def _on_polarization_change(self) -> None:
        """Update the active polarisation channel and trigger a redraw.

        copol: coherent superposition of complex co-pol element patterns.
        cross: coherent superposition of complex cross-pol element patterns.
        total: orthogonal power sum  |AF_copol|² + |AF_xpol|²  (the two
               polarization components are orthogonal so their powers add).
        """
        self._active_polarization = self._polarization_var.get()
        self._recompute_and_redraw()

    def _on_display_mode_change(self) -> None:
        """Apply a change to display mode or dBi range entries.

        Validates the dBi min/max entries (must be floats and max > min),
        flips the colorbar label/ticks to match the chosen mode, and triggers
        a redraw. Invalid entries highlight red and abort the redraw.
        """
        mode = self._display_mode_var.get()

        if mode == DISPLAY_ABSOLUTE:
            # Validate the user-supplied dBi range before applying it.
            try:
                dbi_min = float(self._dbi_min_var.get())
                dbi_max = float(self._dbi_max_var.get())
                if dbi_max <= dbi_min:
                    raise ValueError("max (dBi) must be greater than min (dBi).")
            except ValueError as exc:
                self._dbi_min_entry.config(bg=ENTRY_ERROR_COLOR)
                self._dbi_max_entry.config(bg=ENTRY_ERROR_COLOR)
                self._set_status(f"Invalid dBi range: {exc}")
                return

            self._dbi_min_entry.config(bg=ENTRY_BG_COLOR)
            self._dbi_max_entry.config(bg=ENTRY_BG_COLOR)

            self._cbar.set_label(COLORBAR_LABEL_ABSOLUTE)
            # Span the range with up to 5 evenly spaced ticks.
            self._cbar.set_ticks(np.linspace(dbi_min, dbi_max, 5))
        else:
            self._cbar.set_label(COLORBAR_LABEL_RELATIVE)
            self._cbar.set_ticks(list(range(-DYNAMIC_RANGE_DB, 1, 10)))

        self._recompute_and_redraw()

    # ── CSV PARSING ───────────────────────────────────────────────

    def _parse_weights_csv(self, filepath: str) -> np.ndarray:
        """Parse a weights CSV and return complex weights.

        Compatible with the CSV format written by run_optimization.py.
        Required columns: amplitude (linear), phase_deg (degrees).

        Args:
            filepath (str): Path to the CSV file.

        Returns:
            np.ndarray: Complex weights, shape (N_elements,). Units: dimensionless V/V.

        Raises:
            KeyError: If required column headers are missing.
            ValueError: If the row count does not match the loaded element count.
        """
        amplitudes_linear: list[float] = []
        phases_deg: list[float] = []

        with open(filepath, newline="") as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                if CSV_COL_AMPLITUDE not in row:
                    raise KeyError(
                        f"CSV is missing the '{CSV_COL_AMPLITUDE}' column."
                    )
                if CSV_COL_PHASE_DEG not in row:
                    raise KeyError(
                        f"CSV is missing the '{CSV_COL_PHASE_DEG}' column."
                    )
                amplitudes_linear.append(float(row[CSV_COL_AMPLITUDE]))
                phases_deg.append(float(row[CSV_COL_PHASE_DEG]))

        n_rows = len(amplitudes_linear)
        if n_rows != self._n_elements:
            raise ValueError(
                f"CSV has {n_rows} row(s) but {self._n_elements} element(s) are "
                "loaded. Ensure the CSV matches the element_patterns_dir in config."
            )

        # Convert amplitude + phase (degrees) pairs to complex weights V/V.
        # [MATLAB] use complex(amp .* cosd(phase_deg), amp .* sind(phase_deg))
        phases_rad = np.array(phases_deg) * (np.pi / 180.0)
        amplitudes = np.array(amplitudes_linear)
        weights_complex = amplitudes * np.exp(1j * phases_rad)
        return weights_complex

    # ── RECOMPUTE AND REDRAW ──────────────────────────────────────

    def _recompute_and_redraw(self) -> None:
        """Compute array factor and metrics, then refresh all displays.

        Called on every weight, directive, polarisation, or display-mode
        change. Runs synchronously on the tkinter main thread; numpy
        operations keep latency low for the (N_theta × N_phi) grid.
        """
        directives = self._get_directives_from_ui()

        # ── 1) Compute both copol and cross AFs for CST-convention directivity ──
        # Both stacks are always computed so that directivity can be normalised by
        # total radiated power P_copol + P_cross (IEEE Std 149 partial directivity),
        # matching the denominator CST uses when displaying component patterns.
        # The "active polarisation" controls which power is shown in the heatmap;
        # normalisation always uses the combined power.
        #
        # Power-normalise weights to match the optimizer's cost_fn convention:
        #   w_norm = w / ||w||₂
        # Directivity (|AF|² / P_total) is scale-invariant, so the displayed
        # pattern is unchanged; only Total J becomes consistent with the optimizer.
        # [MATLAB] w_norm = weights_complex / sqrt(max(sum(abs(weights_complex).^2), 1e-30));
        _w_power = float(np.sum(np.abs(self._weights_complex) ** 2))
        _weights_norm = self._weights_complex / np.sqrt(max(_w_power, 1e-30))

        af_copol = compute_array_factor(_weights_norm, self._element_patterns_copol)
        af_cross = compute_array_factor(_weights_norm, self._element_patterns_cross)

        # Spherical integral for total radiated power: P = Σ|AF|² sin(θ) Δθ Δφ
        # [MATLAB] sin_t = sind(theta_deg)'; p = sum(sum((abs(af).^2) .* sin_t)) * dth * dph;
        theta_rad  = np.deg2rad(self._theta_deg)
        dtheta_rad = float(np.diff(theta_rad).mean()) if len(theta_rad) > 1 else np.pi
        dphi_rad   = float(np.deg2rad(np.diff(self._phi_deg).mean())) if len(self._phi_deg) > 1 else 2.0 * np.pi
        sin_theta  = np.sin(theta_rad)

        def _spherical_power(af):
            return float(np.sum(np.abs(af) ** 2 * sin_theta[:, np.newaxis])) * dtheta_rad * dphi_rad

        p_total = _spherical_power(af_copol) + _spherical_power(af_cross)

        # ── 2) Power grid + metrics AF, by polarisation channel ──────
        # copol / cross: |AF|² of a single coherent superposition.
        # total:         |AF_copol|² + |AF_xpol|² (orthogonal channels).
        # For metrics we hand a "pseudo-AF" whose |·|² equals the power grid;
        # phase is irrelevant because evaluate_metrics() only ever uses |AF|².
        if self._active_polarization == POLARIZATION_CROSS:
            power_linear_grid    = np.abs(af_cross) ** 2
            metrics_array_factor = af_cross
        elif self._active_polarization == POLARIZATION_TOTAL:
            # Orthogonal polarizations add in power, not field.
            # [MATLAB] power_total = abs(af_copol).^2 + abs(af_cross).^2;
            power_linear_grid    = np.abs(af_copol) ** 2 + np.abs(af_cross) ** 2
            metrics_array_factor = np.sqrt(power_linear_grid).astype(complex)
        else:  # POLARIZATION_COPOL (default)
            power_linear_grid    = np.abs(af_copol) ** 2
            metrics_array_factor = af_copol

        # ── 3) Build display grid + colorbar limits by display mode ──
        # Always compute the absolute dBi grid — reused for the hover readout
        # regardless of which display mode is active.
        # [MATLAB] dbi_grid = _compute_directivity_dbi_grid(metrics_af, theta, phi, p_total);
        self._dbi_grid = _compute_directivity_dbi_grid(
            metrics_array_factor, self._theta_deg, self._phi_deg,
            normalizer_power=p_total,
        )

        display_mode = self._display_mode_var.get()
        if display_mode == DISPLAY_ABSOLUTE:
            grid_display = self._dbi_grid
            # Parse the user-supplied clim; fall back to defaults on bad input.
            try:
                clim_min = float(self._dbi_min_var.get())
                clim_max = float(self._dbi_max_var.get())
                if clim_max <= clim_min:
                    raise ValueError
            except ValueError:
                clim_min, clim_max = DEFAULT_DBI_MIN, DEFAULT_DBI_MAX
        else:
            # Relative-to-peak dB: normalise so the live peak sits at 0 dB.
            power_db_grid = 10.0 * np.log10(
                np.maximum(power_linear_grid, LOG10_EPSILON)
            )
            grid_display = np.clip(
                power_db_grid - np.max(power_db_grid),
                -DYNAMIC_RANGE_DB, 0.0,
            )
            clim_min, clim_max = -float(DYNAMIC_RANGE_DB), 0.0

        self._update_heatmap(grid_display, clim_min, clim_max)
        self._update_directive_overlays(directives)

        # ── 4) Update live metric labels ─────────────────────────────
        # Pass CST-convention total-power normaliser so that the metric panel
        # shows the same partial directivity that CST displays.
        metrics = evaluate_metrics(
            self._element_patterns_copol,
            self._theta_deg,
            self._phi_deg,
            self._weights_complex,
            directives,
            [],
            precomputed_array_factor=metrics_array_factor,
            normalizer_power=p_total,
        )
        self._update_metrics_labels(metrics)

    def _update_heatmap(self, grid: np.ndarray,
                        clim_min: float, clim_max: float) -> None:
        """Push the display grid to the mesh and update the colorbar limits.

        Args:
            grid (np.ndarray): Display grid, shape (N_theta, N_phi).
                Either normalised-to-peak dB or absolute dBi.
            clim_min (float): Lower colorbar limit.
            clim_max (float): Upper colorbar limit.
        """
        self._last_display_grid = grid
        # QuadMesh.set_array expects a flattened (row-major) 1-D array.
        # [MATLAB] set(mesh_handle, 'CData', grid)
        self._mesh.set_array(grid.ravel())
        self._mesh.set_clim(clim_min, clim_max)
        self._canvas.draw_idle()

    def _update_directive_overlays(self, directives: list[dict]) -> None:
        """Remove old directive artists and redraw for the current directive list.

        Uses the same extended-grid physical masks as the cost function so that
        phi 0°/360° wrap-around and theta pole-crossing are displayed correctly.

        Args:
            directives (list[dict]): Directive dicts with keys type (str),
                theta (float, deg), phi (float, deg), theta_width / phi_width
                or width (float, deg).
        """
        # Remove all previously drawn artists.
        # ContourSet removal API changed in matplotlib 3.8; try .remove() first
        # and fall back to iterating .collections for older versions.
        for artist in self._directive_artists:
            try:
                artist.remove()
            except (ValueError, NotImplementedError, AttributeError):
                if hasattr(artist, "collections"):
                    for coll in artist.collections:
                        try:
                            coll.remove()
                        except (ValueError, NotImplementedError):
                            pass
        self._directive_artists.clear()

        if not directives:
            self._canvas.draw_idle()
            return

        phys_masks = build_directive_physical_masks(
            self._theta_deg, self._phi_deg, directives
        )

        for d, phys_mask in zip(directives, phys_masks):
            color = PEAK_OVERLAY_COLOR if d["type"] == "peak" else NULL_OVERLAY_COLOR

            if phys_mask.any():
                mask_f = phys_mask.astype(float)

                # Semi-transparent fill using a two-stop colormap:
                # 0 (outside window) → fully transparent, 1 (inside) → color + alpha.
                # pcolormesh returns a QuadMesh whose .remove() is reliable across
                # all matplotlib versions (unlike ContourSet in older releases).
                # [MATLAB] use imagesc / patch overlay for mask visualization
                rgba = mcolors.to_rgba(color)
                cmap_ov = mcolors.LinearSegmentedColormap.from_list(
                    "", [(0, 0, 0, 0), (*rgba[:3], WINDOW_ALPHA)]
                )
                mesh_ov = self._ax.pcolormesh(
                    self._phi_deg, self._theta_deg, mask_f,
                    cmap=cmap_ov, vmin=0, vmax=1, zorder=4, shading="auto",
                )
                self._directive_artists.append(mesh_ov)

                # Solid border outline at the mask boundary.
                cs = self._ax.contour(
                    self._phi_deg, self._theta_deg, mask_f,
                    levels=[0.5], colors=[color],
                    linewidths=PATCH_LINEWIDTH, zorder=5,
                )
                # Store contour line collections for removal; handle both old
                # matplotlib (.collections) and new (ContourSet.remove()).
                try:
                    self._directive_artists.extend(cs.collections)
                except AttributeError:
                    self._directive_artists.append(cs)

            # Cross marker at directive target centre.
            # [MATLAB] plot(d_phi, d_theta, '+', 'MarkerSize', 10, 'Color', color);
            scatter = self._ax.scatter(
                d["phi"], d["theta"],
                marker="+", s=SCATTER_MARKER_SIZE,
                color=color, linewidths=PATCH_LINEWIDTH, zorder=6,
            )
            self._directive_artists.append(scatter)

        self._canvas.draw_idle()

    def _update_metrics_labels(self, metrics: dict) -> None:
        """Refresh all fixed metric labels and inline directive result labels.

        Args:
            metrics (dict): Return value of evaluate_metrics().
        """
        total_cost  = metrics.get("total_cost",  None)
        global_peak = metrics.get("global_peak_dbi", None)
        peak_theta  = metrics.get("global_peak_theta_deg", None)
        peak_phi    = metrics.get("global_peak_phi_deg",   None)
        hpbw_th     = metrics.get("hpbw_theta_deg", None)
        hpbw_ph     = metrics.get("hpbw_phi_deg",   None)

        self._label_cost.config(
            text=f"{total_cost:.4e}" if total_cost is not None else "—"
        )
        self._label_peak.config(
            text=f"{global_peak:+.2f} dBi" if global_peak is not None else "—"
        )
        self._label_peak_angle.config(
            text=(f"θ={peak_theta:.1f}° φ={peak_phi:.1f}°"
                  if peak_theta is not None else "—")
        )
        self._label_hpbw.config(
            text=(f"θ:{hpbw_th:.1f}° φ:{hpbw_ph:.1f}°"
                  if hpbw_th is not None else "—")
        )

        # Update inline result labels in the directive table.
        # Walk directive_rows in order, skipping invalid rows (same logic as
        # _get_directives_from_ui), and map the i-th valid row to directive_metrics[i].
        directive_metrics = metrics.get("directive_metrics", [])
        metrics_idx = 0
        for row_vars in self._directive_rows:
            try:
                d_type = row_vars["type_var"].get().strip()
                if d_type not in ("peak", "null"):
                    row_vars["metric_label"].config(text="—", fg="#222222")
                    continue
                float(row_vars["theta"].get())
                float(row_vars["phi"].get())
                float(row_vars["theta_width"].get())
                float(row_vars["phi_width"].get())
                float(row_vars["weight"].get())
            except (ValueError, KeyError):
                row_vars["metric_label"].config(text="—", fg="#222222")
                continue

            if metrics_idx < len(directive_metrics):
                dm = directive_metrics[metrics_idx]
                if dm["type"] == "peak":
                    text  = f"{dm['gain_dbi']:+.2f} dBi"
                    color = PEAK_OVERLAY_COLOR
                else:
                    text  = f"{dm['gain_dbi']:+.2f} ({dm['null_depth_db']:.1f})"
                    color = NULL_OVERLAY_COLOR
                row_vars["metric_label"].config(text=text, fg=color)
                metrics_idx += 1
            else:
                row_vars["metric_label"].config(text="—", fg="#222222")

    # ── MOUSE / POV HANDLERS ──────────────────────────────────────

    def _on_mouse_move(self, event) -> None:
        """Update the in-axes cursor annotation with position and pattern value."""
        if event.inaxes != self._ax or event.xdata is None or self._last_display_grid is None:
            return
        phi_idx   = int(np.argmin(np.abs(self._phi_deg   - event.xdata)))
        theta_idx = int(np.argmin(np.abs(self._theta_deg - event.ydata)))
        val  = self._last_display_grid[theta_idx, phi_idx]
        unit = "dBi" if self._display_mode_var.get() == DISPLAY_ABSOLUTE else "dB"
        self._hover_text.set_text(
            f"θ = {self._theta_deg[theta_idx]:.1f}°   "
            f"φ = {self._phi_deg[phi_idx]:.1f}°   "
            f"D = {val:+.2f} {unit}"
        )
        self._hover_text.set_visible(True)
        self._canvas.draw_idle()

    def _on_mouse_leave(self, _) -> None:
        """Hide the cursor annotation when the mouse leaves the axes."""
        self._hover_text.set_visible(False)
        self._canvas.draw_idle()

    # ── HELPERS ───────────────────────────────────────────────────

    def _get_directives_from_ui(self) -> list[dict]:
        """Read the current directive table state as a list of directive dicts.

        Rows with parse errors (incomplete numeric fields) are silently skipped
        to allow the user to keep typing without breaking the live update.

        Returns:
            list[dict]: Directive dicts for cost_function / evaluate_metrics.
                Keys: type (str), theta (float, deg), phi (float, deg),
                width (float, deg), weight (float, dimensionless).
        """
        directives: list[dict] = []
        for row_vars in self._directive_rows:
            try:
                d_type = row_vars["type_var"].get().strip()
                if d_type not in ("peak", "null"):
                    continue
                directive = {
                    "type":        d_type,
                    "theta":       float(row_vars["theta"].get()),
                    "phi":         float(row_vars["phi"].get()),
                    "theta_width": float(row_vars["theta_width"].get()),
                    "phi_width":   float(row_vars["phi_width"].get()),
                    "weight":      float(row_vars["weight"].get()),
                }
                directives.append(directive)
            except ValueError:
                # Skip rows where a numeric field is still being typed
                pass
        return directives

    def _sync_entries_to_weights(self) -> None:
        """Synchronise all Entry and Slider widgets to match weights_complex.

        Called after Solo or Uniform to keep the displayed values consistent
        with the internal weight state.
        """
        # Guard prevents the slider set() calls below from triggering the slider
        # command callbacks, which would cause a redundant recompute.
        self._syncing_weight_display = True
        for n in range(self._n_elements):
            amplitude_linear = abs(self._weights_complex[n])
            # np.angle returns phase in the range (-π, π]; convert to degrees
            phase_deg = np.angle(self._weights_complex[n], deg=True)

            self._amp_vars[n].set(f"{amplitude_linear:.4f}")
            self._phase_vars[n].set(f"{phase_deg:.2f}")
            self._amp_sliders[n].set(amplitude_linear)
            self._phase_sliders[n].set(phase_deg)

            # Clear any error highlighting from a previous invalid entry
            self._amp_entries[n].config(bg=ENTRY_BG_COLOR)
            self._phase_entries[n].config(bg=ENTRY_BG_COLOR)
        self._syncing_weight_display = False

    def _set_status(self, message: str) -> None:
        """Update the bottom status bar with a new message.

        Args:
            message (str): Text to display.
        """
        self._status_label.config(text=f"  {message}")

    def _update_status_count(self) -> None:
        """Refresh the status bar element-and-directive count."""
        self._set_status(
            f"{self._n_elements} elements loaded · "
            f"{len(self._directive_rows)} directives active"
        )

    # ── RUN ───────────────────────────────────────────────────────

    def run(self) -> None:
        """Start the tkinter event loop (blocks until the window is closed)."""
        self._root.mainloop()


# ────────────────────────── ENTRY POINT ───────────────────────────


def main() -> None:
    """Parse CLI arguments and launch the interactive weight tuner."""
    parser = argparse.ArgumentParser(
        description=(
            "Interactive antenna array weight tuner with live 2-D pattern "
            "and metric feedback."
        )
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        metavar="PATH",
        help="Path to config.yaml (default: %(default)s)",
    )
    args = parser.parse_args()

    tuner = ManualWeightsTuner(config_path=Path(args.config))
    tuner.run()


if __name__ == "__main__":
    main()
