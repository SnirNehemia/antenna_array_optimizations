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
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Rectangle

import yaml

# Add project root to sys.path so src/ modules are importable
# [MATLAB] MATLAB uses addpath(); no sys.path equivalent
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.io.cst_parser import load_element_patterns
from src.cost.cost_function import compute_array_factor
from src.metrics.metrics import evaluate_metrics

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

# Supported polarisation modes (must match keys in element pattern dict)
POLARIZATION_COPOL = "copol"
POLARIZATION_TOTAL = "total"

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

        # Pre-compute stacked element patterns for both supported polarisations.
        # copol: complex co-pol field — E_complex = copol_abs · exp(j · copol_phase)
        # total: total E-field magnitude used as real-valued element patterns;
        #        no phase info for the total field is available in the CST export,
        #        so it is cast to complex with zero imaginary part.
        # [MATLAB] avoid list comprehension; use a for-loop to build a 3-D array
        self._element_patterns_copol: np.ndarray = np.stack(
            [p["E_complex"] for p in element_patterns], axis=0
        )
        self._element_patterns_total: np.ndarray = np.stack(
            [p["E_abs"].astype(complex) for p in element_patterns], axis=0
        )
        # Active pattern stack — switched by the polarisation combobox
        self._element_patterns_stacked: np.ndarray = self._element_patterns_copol

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

        # Matplotlib artists added for directive overlays (cleared on each redraw)
        self._directive_artists: list = []

        # Directive table rows; each is a plain dict keyed by field name → StringVar
        # plus "_frame" → tk.Frame (used for identity-based removal)
        self._directive_rows: list[dict] = []

        # Tkinter widgets created dynamically for per-directive metric lines
        self._per_directive_metric_widgets: list[tk.Widget] = []

        # Frames set during UI build (assigned in _build_* methods)
        self._directives_inner_frame: tk.Frame
        self._metrics_inner_frame: tk.Frame
        self._ax: object            # matplotlib Axes
        self._mesh: object          # QuadMesh returned by pcolormesh
        self._canvas: FigureCanvasTkAgg
        self._status_label: tk.Label
        self._label_cost: tk.Label
        self._label_peak: tk.Label
        self._label_pnr: tk.Label

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
        self._build_status_bar()

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
            values=[POLARIZATION_COPOL, POLARIZATION_TOTAL],
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

        # 2-D heatmap: azimuth φ on x-axis, elevation θ on y-axis
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

        cbar = fig.colorbar(
            self._mesh, ax=self._ax, label="dB (normalised to peak)"
        )
        cbar_ticks = list(range(-DYNAMIC_RANGE_DB, 1, 10))
        cbar.set_ticks(cbar_ticks)

        self._ax.set_xlabel("Azimuth  φ (°)")
        self._ax.set_ylabel("Elevation  θ (°)")
        self._ax.set_title("Array Factor — manual weights", fontsize=10)
        self._ax.grid(
            True, alpha=GRID_ALPHA, color="white", linewidth=GRID_LINEWIDTH
        )

        # Fix axes to full-sphere coverage regardless of data resolution.
        self._ax.set_xlim(0.0, 360.0)
        self._ax.set_ylim(0.0, 180.0)

        fig.tight_layout()

        self._canvas = FigureCanvasTkAgg(fig, master=lf)
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

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
            ("Type", 6), ("θ (°)", 6), ("φ (°)", 6), ("W (°)", 6), ("wt", 5)
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
        self._label_cost = self._make_metric_label("Total J:", "—")
        self._label_peak = self._make_metric_label("Global peak:", "—")
        self._label_pnr  = self._make_metric_label("Peak-to-null:", "—")

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

    def _build_status_bar(self) -> None:
        """Build the status bar along the bottom of the window."""
        self._status_label = tk.Label(
            self._root,
            text="",
            font=LABEL_FONT,
            bg=PANEL_BG,
            fg="#555555",
            anchor="w",
            relief=tk.SUNKEN,
            bd=1,
        )
        self._status_label.pack(side=tk.BOTTOM, fill=tk.X)

    # ── DIRECTIVE TABLE MANAGEMENT ────────────────────────────────

    def _add_directive_row(self, directive_dict: dict | None) -> None:
        """Append a new row to the directives table and trigger a redraw.

        Args:
            directive_dict (dict | None): Pre-populated values from config, or
                None to use project defaults. Expected keys: type, theta (deg),
                phi (deg), width (deg), weight.
        """
        # Fall back to default values for any missing field
        defaults = {
            "type":   DEFAULT_DIRECTIVE_TYPE,
            "theta":  DEFAULT_DIRECTIVE_THETA_DEG,
            "phi":    DEFAULT_DIRECTIVE_PHI_DEG,
            "width":  DEFAULT_DIRECTIVE_WIDTH_DEG,
            "weight": DEFAULT_DIRECTIVE_WEIGHT,
        }
        d = directive_dict if directive_dict is not None else {}

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

        # Numeric fields: theta, phi, width, weight (all in degrees or dimensionless)
        numeric_field_specs = [
            ("theta",  d.get("theta",  defaults["theta"]),  6),
            ("phi",    d.get("phi",    defaults["phi"]),    6),
            ("width",  d.get("width",  defaults["width"]),  6),
            ("weight", d.get("weight", defaults["weight"]), 5),
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
        """Switch the active element pattern stack based on the polarisation selector.

        copol: complex co-pol patterns (E_complex from parser).
        total: real total-field magnitude patterns (E_abs cast to complex;
               no cross-pol phase info is available in the CST export).
        """
        if self._polarization_var.get() == POLARIZATION_COPOL:
            self._element_patterns_stacked = self._element_patterns_copol
        else:
            # Total field: magnitude-only real-valued patterns cast to complex
            self._element_patterns_stacked = self._element_patterns_total
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

        Called on every weight or directive change. Runs synchronously on
        the tkinter main thread; numpy operations keep latency low.
        """
        directives = self._get_directives_from_ui()

        # Array factor: coherent superposition of weighted element patterns.
        # AF(θ,φ) = Σ_n  w_n · E_n(θ,φ),  shape (N_theta, N_phi), V/m.
        array_factor = compute_array_factor(
            self._weights_complex, self._element_patterns_stacked
        )

        # Power density in linear scale (V/m)².
        # Project convention: all optimisation math uses linear magnitude, not dB.
        power_linear_grid = np.abs(array_factor) ** 2

        # Convert to dB, normalise so the peak sits at 0 dB, clip to dynamic range.
        power_db_grid = 10.0 * np.log10(
            np.maximum(power_linear_grid, LOG10_EPSILON)
        )
        peak_db = np.max(power_db_grid)
        grid_norm_db = np.clip(power_db_grid - peak_db, -DYNAMIC_RANGE_DB, 0.0)

        self._update_heatmap(grid_norm_db)
        self._update_directive_overlays(directives)

        # Evaluate directivity-based metrics; pass empty cost_history (not applicable
        # in manual mode; only used for iteration_count which is not displayed here)
        metrics = evaluate_metrics(
            self._element_patterns_stacked,
            self._theta_deg,
            self._phi_deg,
            self._weights_complex,
            directives,
            [],
        )
        self._update_metrics_labels(metrics)

    def _update_heatmap(self, grid_norm_db: np.ndarray) -> None:
        """Push new data into the existing pcolormesh without rebuilding the axes.

        Args:
            grid_norm_db (np.ndarray): Normalised power in dB, shape (N_theta, N_phi).
                Values in the range [-DYNAMIC_RANGE_DB, 0]. Units: dB.
        """
        # QuadMesh.set_array expects a flattened (row-major) 1-D array.
        # [MATLAB] set(mesh_handle, 'CData', grid_norm_db)
        self._mesh.set_array(grid_norm_db.ravel())
        self._canvas.draw_idle()

    def _update_directive_overlays(self, directives: list[dict]) -> None:
        """Remove old directive artists and redraw for the current directive list.

        Args:
            directives (list[dict]): Directive dicts with keys type (str),
                theta (float, deg), phi (float, deg), width (float, deg).
        """
        # Remove all previously drawn patches and scatter markers
        for artist in self._directive_artists:
            try:
                artist.remove()
            except (ValueError, NotImplementedError):
                pass
        self._directive_artists.clear()

        for d in directives:
            color = (
                PEAK_OVERLAY_COLOR if d["type"] == "peak" else NULL_OVERLAY_COLOR
            )
            # Half angular width of the enforcement / suppression window, deg
            half_width_deg = d["width"] / 2.0

            # Rectangle marks the full angular window on the θ-φ heatmap
            rect = Rectangle(
                (d["phi"] - half_width_deg, d["theta"] - half_width_deg),
                2.0 * half_width_deg,
                2.0 * half_width_deg,
                linewidth=PATCH_LINEWIDTH,
                edgecolor=color,
                facecolor=color,
                alpha=WINDOW_ALPHA,
                zorder=4,
            )
            self._ax.add_patch(rect)
            self._directive_artists.append(rect)

            # Cross marker at exact target (φ, θ) coordinate
            scatter = self._ax.scatter(
                d["phi"],
                d["theta"],
                marker="+",
                s=SCATTER_MARKER_SIZE,
                color=color,
                linewidths=PATCH_LINEWIDTH,
                zorder=5,
            )
            self._directive_artists.append(scatter)

        self._canvas.draw_idle()

    def _update_metrics_labels(self, metrics: dict) -> None:
        """Refresh all fixed and per-directive metric labels.

        Args:
            metrics (dict): Return value of evaluate_metrics().
                Expected keys: total_cost, global_peak_dbi,
                peak_to_null_ratio_db, directive_metrics.
        """
        total_cost    = metrics.get("total_cost", None)
        global_peak   = metrics.get("global_peak_dbi", None)
        pnr_db        = metrics.get("peak_to_null_ratio_db", None)

        self._label_cost.config(
            text=f"{total_cost:.4e}" if total_cost is not None else "—"
        )
        self._label_peak.config(
            text=f"{global_peak:+.2f} dBi" if global_peak is not None else "—"
        )
        self._label_pnr.config(
            text=f"{pnr_db:.2f} dB" if pnr_db is not None else "—"
        )

        # Destroy existing per-directive rows and rebuild from fresh metrics
        for widget in self._per_directive_metric_widgets:
            widget.destroy()
        self._per_directive_metric_widgets.clear()

        directive_metrics = metrics.get("directive_metrics", [])
        if directive_metrics:
            sep = ttk.Separator(self._metrics_inner_frame, orient="horizontal")
            sep.pack(fill=tk.X, pady=3)
            self._per_directive_metric_widgets.append(sep)

        for dm in directive_metrics:
            gain_dbi  = dm.get("gain_dbi", None)
            d_type    = dm.get("type", "?").capitalize()
            d_theta   = dm.get("theta_deg", 0.0)
            d_phi     = dm.get("phi_deg", 0.0)

            row = tk.Frame(self._metrics_inner_frame, bg=BG_COLOR)
            row.pack(fill=tk.X, pady=1)

            tk.Label(
                row,
                text=f"{d_type} θ={d_theta:.1f}° φ={d_phi:.1f}°:",
                font=LABEL_FONT,
                bg=BG_COLOR,
                anchor="w",
                width=22,
            ).pack(side=tk.LEFT)

            value_text = f"{gain_dbi:+.2f} dBi" if gain_dbi is not None else "—"
            tk.Label(
                row,
                text=value_text,
                font=METRIC_FONT,
                bg=BG_COLOR,
                anchor="w",
                fg="#222222",
            ).pack(side=tk.LEFT)

            self._per_directive_metric_widgets.append(row)

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
                    "type":   d_type,
                    "theta":  float(row_vars["theta"].get()),
                    "phi":    float(row_vars["phi"].get()),
                    "width":  float(row_vars["width"].get()),
                    "weight": float(row_vars["weight"].get()),
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
