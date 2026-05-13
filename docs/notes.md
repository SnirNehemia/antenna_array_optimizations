# Development Notes

## Physics & Modelling Assumptions

### Far-field assumption
All element patterns exported from CST are far-field patterns. The array model
therefore uses the far-field superposition principle (no near-field coupling
between elements). Mutual coupling effects are captured implicitly only if CST
was run with the full array excited (i.e., active element patterns). If patterns
are exported from single-element simulations in isolation, mutual coupling is
**not** accounted for. Clarify this when loading patterns from a new dataset.

### Co-polarization field as the optimization target
The complex element pattern used in the array factor is constructed from the
**co-polarization** component (columns `Abs(Copol)` + `Phase(Copol)`). This is
appropriate when the array is intended to be used with a fixed polarization.
Cross-polarization optimization is a future extension; a config flag
(`use_total_field: bool`) can switch to `Abs(E)` if total radiated power matters
more than polarization purity.

### Array factor model (no steering vector)
The optimizer works directly with sampled complex element patterns on a discrete
angular grid — it does **not** assume a uniform linear array (ULA) model or an
analytical steering vector `exp(j k d sin θ)`. This makes it general: it handles
arbitrary element positions, non-isotropic patterns, and irregular arrays.
The trade-off is that gradient computation scales with the number of angular
samples × number of elements (rather than being analytically derived).

### 2D vs 3D patterns
The current pipeline stores full 3D patterns (Theta × Phi grid) but the initial
optimization runs target a **2D azimuth cut** (fixed Phi, varying Theta). This
is the simplest starting point. Extension to full 3D optimization is possible
by expanding the angular window of each directive to cover all Phi slices.

### Angular window for directives
Each directive specifies an angular `width` in degrees. The cost function
accumulates the squared field magnitude over **all grid points within ±width/2**
of the target angle. This soft-window approach avoids discontinuities in the
gradient. For very fine angular targets (width < angular resolution), the
nearest grid point is used.

---

## Numerical Conventions

### Optimization variable layout
The optimizer uses a real-valued vector `x` of length `2N`:
```
x = [Re(w_0), Im(w_0), Re(w_1), Im(w_1), ..., Re(w_{N-1}), Im(w_{N-1})]
```
This is the standard way to feed complex-valued problems to real-valued solvers
(L-BFGS-B, BFGS, Nelder-Mead). MATLAB's `fminunc`/`fmincon` use the same trick.

### Phase representation
- Internal computation: **radians** (numpy default)
- Display / output CSV: **degrees**
- CST import: Phase columns are in **degrees** → convert on import with `* np.pi / 180`

### dB conversion
```
pattern_dB = 20 * log10(|AF(theta, phi)| / |AF_max|)
```
Use 20·log10 (field quantity), not 10·log10 (power). Normalization to the peak
of the optimized pattern unless absolute gain is requested.

### Cost function sign
L-BFGS-B **minimizes**. Therefore:
- Null directives: positive cost (minimize = suppress)
- Peak directives: **negative** of the desired gain (minimize = maximize)

---

## Known Issues / Open Questions

1. **Local minima**: Gradient descent is not globally optimal. Multi-start
   (random initial weights, run K times, keep best) is implemented but
   not yet validated. Flag if results look unexpected.

2. **Angular grid interpolation**: If a directive targets an angle between
   two grid points, the current implementation uses the nearest grid point.
   Bilinear interpolation on the (Theta, Phi) grid would be more accurate
   but adds complexity — defer until needed.

3. **Element count not yet fixed**: The pipeline is written for N elements but
   the actual array geometry and element count will be defined when the user
   provides all N element pattern files. Until then, the parser is tested on
   a single element file.

4. **Phase(Copol) wraparound**: CST outputs phase in degrees over [0°, 360°]
   or [-180°, 180°] depending on settings. The parser should handle both.
   Currently assumes [0°, 360°]. TODO: add phase unwrapping if needed.

5. **`Ax.Ratio` column**: Currently parsed but not used. Reserved for future
   polarization-diversity optimization.

---

## MATLAB Porting Notes

> Track all Python-specific constructs here so the migration is friction-free.

### Python constructs that need adaptation

| Python construct | MATLAB equivalent | Notes |
|---|---|---|
| `numpy.ndarray` | `double` array (matrix) | MATLAB is 1-indexed; adjust all slice indices |
| `scipy.optimize.minimize(method='L-BFGS-B')` | `fminunc` (unconstrained) or `fmincon` (with bounds) | `fmincon` for amplitude bounds |
| `numpy.exp(1j * x)` | `exp(1i * x)` | MATLAB uses `1i` or `1j` for imaginary unit |
| `numpy.reshape(arr, (M, N))` | `reshape(arr, [M, N])` | Same semantics, different syntax |
| `numpy.linalg.norm` | `norm` | Built-in in MATLAB |
| `f-strings` / `str.format` | `sprintf` | For output formatting |
| `dict` for structured data | `struct` | Replace Python dicts with MATLAB structs |
| List of dicts (directives) | Array of structs | `directives(k).type`, `.angle`, etc. |
| `enumerate` | `for k = 1:N` | MATLAB for-loop over index |
| `np.pi` | `pi` | MATLAB built-in |
| `complex(re, im)` | `re + 1i*im` | Explicit in MATLAB |
| `os.path`, `glob` | `fullfile`, `dir` | Filesystem operations |
| `pytest` | No direct equivalent; use MATLAB's `assert` in scripts | |

### Functions that map 1-to-1

```
cst_parser.py         →  cst_parser.m
cost_function.py      →  cost_function.m
optimizer.py          →  optimizer.m          (wraps fminunc/fmincon)
plotter.py            →  plotter.m            (matplotlib → MATLAB figure/polarplot)
metrics.py            →  metrics.m
run_optimization.py   →  run_optimization.m
```

### Key indexing difference
Python: `x[0]` is the first element.
MATLAB: `x(1)` is the first element.
The weight vector layout `[Re(w_0), Im(w_0), ...]` maps to
`[Re_w(1), Im_w(1), ...]` in MATLAB — adjust all index arithmetic by +1.

### gradient / Jacobian
L-BFGS-B uses finite-difference gradients by default (no need to derive Jacobian
analytically). `fminunc` also supports finite differences via `'FiniteDifferenceType'`
option. For speed, an analytic gradient can be added later in both languages.

---

## Future Work (Out of Scope for Current Phase)

- **Dynamic jammer environment**: Adaptive null steering in response to moving
  interference sources. This will require a real-time update loop and possibly
  a different solver (e.g., convex optimization via CVXPY / MATLAB CVX).
  Keep the cost function and array model modular so they can be called
  incrementally.

- **Full 3D pattern optimization**: Extend directives to specify (theta, phi)
  pairs instead of 1D angle cuts.

- **Cross-polarization suppression**: Add a cross-pol directive type that
  penalizes `Abs(Cross)` in specified angular regions.

- **Mutual coupling correction**: If single-element patterns are used, add a
  Z-matrix or S-matrix correction layer before the array factor computation.

---

## Config Schema Reference

The user-facing configuration is `config.yaml` in the project root.
Key sections and their types:

```yaml
element_patterns_dir: str          # path to folder of CST .txt files
polarization: "copol" | "total"    # which field component to optimize

directives:                        # list of beam-shaping directives
  - type:   "peak" | "null"
    theta:  float                  # degrees, elevation
    phi:    float                  # degrees, azimuth [optional, default 0.0]
    width:  float                  # degrees, angular window
    level:  float | null           # dBi (peak) or dB re peak (null) [optional]
    weight: float                  # cost function lambda_k [optional, default 1.0]

optimizer:
  max_iterations:    int           # L-BFGS-B max iterations
  cost_tolerance:    float         # convergence threshold on J
  n_restarts:        int           # number of random multi-starts
  amplitude_bounds:  [float, float] | null   # [min, max] per element
  phase_only:        bool          # fix amplitudes=1, optimize phase only

output:
  results_dir:            str      # base folder for timestamped output
  plot_cut_type:          "theta_cut" | "phi_cut"
  plot_phi_deg:           float
  plot_theta_deg:         float
  save_polar_plot:        bool
  save_cartesian_plot:    bool
  save_weight_plots:      bool
  save_cost_history_plot: bool
```

If any required key is missing, the code must raise a descriptive `KeyError`
or `ValueError` — never silently fall back to a hardcoded default.

---

## Session Log

> Claude Code must append an entry here at the end of every working session.
> Format shown below. Newest entry at the top.

### Template
```
### YYYY-MM-DD — <one-line summary>
**Implemented**: ...
**Decisions made**: ...
**Open questions / known issues**: ...
```

### 2026-05-13 — Interactive manual weight tuner (manual_weights.py)

**Implemented**:
- `scripts/manual_weights.py` — new standalone interactive GUI. Loads element patterns
  and config directives, displays a live 2-D radiation pattern heatmap (jet colormap,
  normalised to peak, −40 dB dynamic range), and updates it in real time whenever the
  user changes weights or directive targets. Launched with
  `python scripts/manual_weights.py --config config.yaml`.
- UI layout: left panel (~60 % of width) = matplotlib `pcolormesh` canvas embedded via
  `FigureCanvasTkAgg`; right panel = scrollable element-weight controls + directives table
  + live metrics.
- **Element weights panel**: one row per element with a narrow amplitude entry + horizontal
  `ttk.Scale` slider ([0, 2]) and a phase entry + horizontal slider ([−180°, 180°]).
  Entries and sliders are bidirectionally synced via a `_syncing_weight_display` guard flag
  that suppresses recursive callbacks during programmatic updates. Entries commit on
  `<Return>` or `<FocusOut>`; sliders fire immediately on drag.
- **Solo button** (per element): zeros all other elements and sets the selected element to
  1+0j — useful for inspecting individual element patterns.
- **Uniform Weights** button: resets all elements to amplitude = 1.0, phase = 0.0°.
- **Load Weights CSV** button: opens a file dialog and reads the standard
  `weights.csv` format written by `run_optimization.py` (columns: amplitude, phase_deg).
  Row count is validated against the loaded element count.
- **Directives table**: pre-populated from `config["directives"]` on launch. Each row has
  a type combobox (peak/null), θ / φ / width / weight entries, and a × remove button.
  "+ Add Directive" appends a default row. Active directives are overlaid on the heatmap as
  green (peak) or red (null) Rectangle patches + cross markers, matching the style of the
  existing `pattern_2d.png` output.
- **Metrics panel**: displays Total J, Global peak (dBi), Peak-to-null ratio (dB), and one
  gain (dBi) line per active directive. Rebuilt dynamically when directive rows are added or
  removed. All values computed via `evaluate_metrics()` (called with `cost_history=[]`).
- **Polarisation selector** (`Combobox` in toolbar): switches between `copol` (uses
  `E_complex` — co-pol magnitude + phase) and `total` (uses `E_abs` cast to complex,
  magnitude-only real-valued element patterns). Both stacks are pre-computed at load time.
  Note: no phase information is available for the total field in the CST export; the `total`
  mode effectively treats element patterns as in-phase magnitude envelopes.
- **Axis limits** fixed explicitly to θ ∈ [0°, 180°] and φ ∈ [0°, 360°] in both
  `manual_weights.py` (`_build_pattern_panel`) and `src/plot/plotter.py`
  (`save_2d_projection_plot`), so the full-sphere domain is always visible regardless of
  the angular resolution of the loaded data.
- Config file opened with `encoding="utf-8"` (same fix as `run_optimization.py`) to avoid
  `UnicodeDecodeError` on Windows with the default cp1252 codec.

**Decisions made**:
- Heatmap updated in-place via `QuadMesh.set_array(grid.ravel())` + `canvas.draw_idle()`
  rather than rebuilding the axes on each change. This keeps redraw latency low enough for
  slider drag to feel responsive on the 181×360 (16-element) dataset.
- `total` polarisation uses `E_abs.astype(complex)` (zero imaginary part). The alternative
  of computing a coherent total-field complex pattern from copol + cross-pol would require
  the parser to also extract `cross_complex`, which is out of scope for this session.
  Documented as a known limitation in the UI label.
- Slider range for amplitude is [0, 2] (not locked to `amplitude_bounds` from config) to
  avoid a hard dependency on the optimizer config section, which is irrelevant for this tool.
- Directive rows use `<Return>` / `<FocusOut>` bindings (not StringVar trace) to avoid
  recomputing on every keystroke while the user is still typing.

**Open questions / known issues**:
- `total` polarisation mode is an approximation (magnitude-only, no cross-pol phase). A
  rigorous implementation would require `cst_parser.py` to return `cross_complex =
  cross_abs * exp(j * cross_phase_rad)` and a coherent total-field combination. Filed as
  future work.
- The directive width in the overlay Rectangle uses the same `±width/2` box convention as
  the cost function. If the cost function is later changed to a circular window, the overlay
  would need updating accordingly.

### 2026-05-12 — Multi-start convergence plot, single-element inits, power normalization

**Implemented**:
- `cost_function.py` — `cost_fn` now normalizes weights by `||w||₂` before computing the array factor: `AF_norm = AF / ||w||₂`. Physical model: fixed total source power split across elements (power-splitter). The cost is now scale-invariant — inflating all weights by the same factor yields no gain, so the optimizer no longer saturates to the Re/Im box corners (`sqrt(2)` amplitude). Only relative phases and amplitudes matter.
- `optimizer.py` — new `_single_element_initial_x(n_elements, element_idx, mode)`: weight vector with element `element_idx` = 1+0j and all others zero. Used as additional starting points.
- `optimizer.py` — `run_optimizer` now runs two phases: phase 1 = `n_restarts` user-configured runs (uniform + random); phase 2 = one run per element with only that element active. Tracks `all_cost_histories`, `all_run_labels`, `best_run_index` alongside the existing `cost_history` / `result` / `weights_complex`.
- `optimizer.py` — two new optional config keys: `use_uniform_init` (default True, controls whether run 0 uses uniform weights) and `use_single_element_init` (default True, controls whether per-element runs are added). `ValueError` raised if both paths are disabled and no runs would execute.
- `plotter.py` — `save_cost_history_plot` rewritten: accepts `all_cost_histories`, `best_run_index`, `all_run_labels`. All non-best runs drawn as thin grey lines with a single "Other runs (N total)" legend entry; best run drawn bold blue with its init label (e.g. "Best — Element 3 init").
- `save_all_plots` / `run_optimization.py` API updated to pass `all_cost_histories`, `best_run_index`, `all_run_labels` through the call chain. `evaluate_metrics` and `run_report` still use `cost_history` (best run only).
- Minor: `_print_summary` updated to include phi in per-directive printed lines (matching `run_report.txt`). Removed the `# run: <timestamp>` comment header from `weights.csv`.
- `config.yaml` — documented `use_uniform_init` and `use_single_element_init` flags.

**Decisions made**:
- Power normalization is applied inside `cost_fn` (not as a post-processing step), so the optimizer objective is always physically meaningful regardless of initial scale.
- `metrics.py` and `plotter.py` are unaffected: directivity already normalizes by integrating `|AF|²` over the sphere (scale-invariant), so dBi values are unchanged.
- Single-element inits are skipped for `phase_only` mode because all amplitudes are normalized to 1 inside the cost function — the starting amplitude is irrelevant.
- Guard `if best_result is None` added to give a clear error if the user disables all run paths instead of a cryptic `AttributeError`.

**Open questions / known issues**:
- With `amplitude_bounds: [0, 1]` in standard mode, the actual amplitude constraint is still `|Re(w_n)| ≤ 1` and `|Im(w_n)| ≤ 1` (a square, not a disk). Power normalization removes the incentive to saturate to `sqrt(2)`, but the feasible set geometry is still a box. A polar parameterization (`amplitude ∈ [0,1]`, `phase ∈ [0, 2π]`) would be the rigorous fix for a true `|w_n| ≤ 1` constraint.

### 2026-05-12 — Results archive, visualization improvements, cross-section filtering

**Implemented**:
- `run_optimization.py` — new `_save_run_report()`: writes `run_report.txt` to each timestamped results folder with run timestamp, optimizer settings, directives, and per-directive dBi results. Enables side-by-side comparison of two runs.
- `run_optimization.py` — `shutil.copy(args.config, output_dir / "config.yaml")` archives the exact config used in each run.
- `run_optimization.py` — `_save_weights_csv`: prepends a `# run: <timestamp>` comment line to `weights.csv` so each file is self-identifying when browsing multiple result folders.
- `plotter.py` — `save_2d_projection_plot`: replaced single-point scatter markers with `matplotlib.patches.Rectangle` patches showing the full angular window extent (2·width × 2·width), semi-transparent filled + solid border + center cross. Clearly communicates both the target location and window size.
- `plotter.py` — new `_directive_on_cut(directive, output_config)`: replaces `_directive_on_front_half`. Returns `(visible, on_front_half)`. A directive is visible only if its target phi (theta_cut) or target theta (phi_cut) is within `directive["width"]/2` of the cut's fixed angle. Directives that miss the cut are silently skipped in polar and Cartesian plots.
- `plotter.py` — new `_resolve_output_config(output_config, directives)`: resolves `plot_cut_type: "auto"` to `theta_cut` at the first peak directive's phi. Called once at the start of `save_all_plots`; individual plot functions receive an already-resolved config.
- `config.yaml` — `plot_cut_type` changed to `"auto"` (self-setting default). Optimizer stop-criteria comments expanded with the scipy `ftol` formula and guidance on `cost_tolerance` and `n_restarts` usage.

**Decisions made**:
- Directive visibility threshold uses the directive's own half-width (not a fixed angle like 45°), so a narrow directive only appears on a cut that genuinely intersects its window.
- `_resolve_output_config` returns a shallow copy; the caller's dict is never mutated. Individual plot functions remain callable standalone (they just expect a resolved config — "auto" would fall through to phi_cut default, not a crash).
- Run report timestamp is captured at the time of writing (not at start of the run), which is close enough and avoids threading the start time through more function signatures.

**Open questions / known issues**:
- Weights "same always" across runs: confirmed by design — `seed=0` makes multi-start reproducible for a given config. User should vary the config between runs and compare `run_report.txt` files.
- If the user can reproduce a within-run CSV ≠ PNG discrepancy, code inspection shows the same `weights_complex` object is passed to both; the cause would need to be reproduced before a fix can be targeted.

### 2026-05-12 — dBi metrics: directivity-based gain reporting

**Implemented**:
- `src/metrics/metrics.py` — new `_compute_directivity_dbi_grid()`: numerically integrates `|AF|² · sin(θ) · Δθ · Δφ` over the full sphere (rectangle rule). Absolute V/m scale from CST patterns cancels in `4π|AF|²/P_total`, yielding true directivity in dBi.
- `evaluate_metrics` updated: `global_peak_dbi`, per-directive `gain_dbi`, `null_depth_db` now all dBi-based. Removed dead `_linear_to_db()`.
- `run_optimization.py` — `_print_summary` updated to print dBi labels; null lines show both `{gain_dbi:.2f} dBi` and `(depth = {null_depth_db:.2f} dB)`.

**Decisions made**:
- Expected `global_peak_dbi` for a 16-element broadside array is ~14–20 dBi (vs the ~31 dB that appeared with raw `10·log10(|AF|²)` from absolute-V/m patterns).
- Grid spacings Δθ, Δφ derived from `np.diff(...).mean()` to handle non-uniform grids.

**Open questions / known issues**:
- None new.

### 2026-05-11 — Cosmetic plot fixes and 2D projection view

**Implemented**:
- `run_optimization.py`: wall-clock timer around `run_optimizer()`; elapsed time printed in the Results block using `time.perf_counter()`.
- `plotter.py` — `save_polar_plot`: for `theta_cut`, stitches the phi front-half (0→180°) and phi+180° back-half into a full -180°→+180° sweep so the polar plot fills a complete circle. Directive windows placed on whichever phi-half they belong to.
- `plotter.py` — `save_cartesian_plot`: same front/back extraction; x-axis locked to ±180°; normalized to peak = 0 dB (previously showed raw absolute dB). X-axis label notes which phi each side represents.
- Both pattern plots: y/radial axis locked to `[-plot_dynamic_range_db, 0]` dB (shared, configurable, default 40 dB).
- `plotter.py` — new `save_2d_projection_plot()`: `pcolormesh` of the full (N_theta × N_phi) grid normalized to 0 dB, `jet` colormap. Green `*` / red `x` scatter markers at each directive's (phi, theta). Controlled by `save_2d_projection_plot` config flag.
- `config.yaml`: added `plot_dynamic_range_db: 40` and `save_2d_projection_plot: true`.

**Decisions made**:
- Back half of theta_cut uses the actual phi+180° column from the measured data, not a mirror of the front half, so the plot is physically correct for non-symmetric patterns.
- `DEFAULT_DYNAMIC_RANGE_DB = 40` constant added as module-level fallback; config key overrides it at runtime.
- Config schema reference table in this file does not yet reflect the new `plot_dynamic_range_db` and `save_2d_projection_plot` keys — update when schema section is next revisited.

**Open questions / known issues**:
- None new.

### 2026-05-11 — Stages 5–7: metrics, visualization, and full pipeline entry point

**Implemented**:
- `src/metrics/__init__.py` + `src/metrics/metrics.py` (Stage 5) — `evaluate_metrics()` computes per-directive `gain_db` at the nearest grid point, `null_depth_db` relative to the global pattern peak, per-directive cost terms, `peak_to_null_ratio_db`, and `global_peak_db`. Imports only `compute_array_factor` and `angular_window_mask` from `src.cost.cost_function`.
- `src/plot/__init__.py` + `src/plot/plotter.py` (Stage 6) — `save_polar_plot`, `save_cartesian_plot`, `save_weight_plots`, `save_cost_history_plot`, and the `save_all_plots` orchestrator. Headless-safe via `matplotlib.use("Agg")`. Directive windows shaded green (peak) / red (null) in both plot types. Polar pattern normalized to 0 dB at peak.
- `scripts/run_optimization.py` (Stage 7) — full pipeline entry point: reads `config.yaml`, runs Stages 1→4→5→6, writes `weights.csv` (element_index, amplitude, phase_deg, real, imag) and `metrics.json` to a timestamped folder under `results/`. Accepts `--config` flag. Validated end-to-end on 16-element data: converged in 28 iterations, peak 18.89 dB, null depth −49.43 dB, peak-to-null 37.27 dB.

**Decisions made**:
- Metrics reports power at the **exact target grid point** (nearest-index), not the window average used by the cost function. Window average is an optimizer quantity; point value is the physical metric.
- Null depth always referenced to `global_peak_db` (not to a specific peak directive) so it is well-defined even when no peak directive exists.
- Weight plots produced as a **single PNG** with two vertically-stacked subplots (amplitude + phase), controlled by one `save_weight_plots` flag.
- Plotter receives a pre-computed `array_factor_db_grid` from the entry point; it does not import from `src.cost` or `src.optimize`. Module boundary is preserved.
- Polarization guard in entry point: raises `ValueError` with a clear message if `polarization` is anything other than `"copol"`. Prevents silent incorrect results until total-field support is added.
- `sys.path.insert(0, project_root)` added to `scripts/run_optimization.py` so it can be invoked directly as `python scripts/run_optimization.py` without `PYTHONPATH` configuration.
- Config YAML opened with `encoding="utf-8"` to handle box-drawing characters in comments on Windows (default cp1252 codec fails).

**Open questions / known issues**:
- `notes.md` §Numerical Conventions states `20·log10(|AF| / |AF_max|)` for dB conversion. The implementation uses `10·log10(|AF|²)`, which is mathematically identical (`10·log10(|AF|²) = 20·log10(|AF|)`), but normalization differs: code uses absolute power (not peak-normalized) except in the polar plot. No action needed, but the note wording could be clarified.
- Known issue #3 ("element count not yet fixed") is now resolved: parser validated on all 16 elements of `data/Env_1_1/`.
- Config schema reference in notes does not document the `amplitude_only` optimizer key added in a prior session. Should be updated.

### 2026-05-11 — Stages 1–4: CST parser, cost function, and L-BFGS-B optimizer

**Implemented**:
- `src/__init__.py`, `src/io/__init__.py`, `src/cost/__init__.py`, `src/optimize/__init__.py` — package init files.
- `src/io/cst_parser.py` (Stage 1) — `_detect_grid_shape` (auto-detects Theta/Phi resolution via `np.unique`), `_extract_element_index` (regex on filename `\[(\d+)\]`), `parse_cst_file` (skips 2-line header, parses 8 columns, reshapes flat data to 2D grid, builds `E_complex` from Copol magnitude + phase), `load_element_patterns` (loads directory, sorts by element index, validates shared grid shape). Validated on 16 real CST files in `data/Env_1_1/` (181×360 grid, 1° resolution).
- `src/cost/cost_function.py` (Stage 3) — `x_to_weights` / `weights_to_x` (Re/Im interleaved 2N encoding), `compute_array_factor` (coherent superposition), `angular_window_mask` (boolean 2D mask, full-width interpretation), `_directive_cost` (peak: negative mean power; null: positive mean power), `build_cost_function` (pre-computes masks, returns closure; handles `"standard"`, `"phase_only"`, `"amplitude_only"` modes).
- `src/optimize/optimizer.py` (Stage 4) — `_build_lbfgsb_bounds`, `_uniform_initial_x`, `_random_initial_x`, `_run_single_optimization` (callback records cost history), `run_optimizer` (multi-start L-BFGS-B, seeded `np.random.default_rng(seed=0)`, returns best by `result.fun`).
- `config.yaml` — added `amplitude_only: false` to the `optimizer:` section alongside the existing `phase_only` flag.

**Decisions made**:
- **Reshape direction (theta-fast axis)**: `STYLE.md` contained a contradictory example (`flat.reshape(n_theta, n_phi)` with comment "Theta varies fastest"). Actual data rows were read directly to confirm theta is the fast (inner) loop. Correct reshape: `flat.reshape(n_phi, n_theta).T` → shape `(N_theta, N_phi)`.
- **amplitude_only mode**: Added at user request after the initial plan only covered `phase_only`. Uses a length-N real variable vector where `x[n]` is the amplitude of element n directly (phase = 0). Bounds `[a_min, a_max]` are exact (no Re/Im approximation). Decoded via `x.astype(complex)`.
- **Variable encoding lives in `cost_function.py`**: `x_to_weights` and `weights_to_x` are defined there because the cost function owns the encoding contract; the optimizer imports them rather than duplicating the logic.
- **Masks pre-computed in closure**: Angular window masks are constant across all optimizer iterations; computing them once inside `build_cost_function` before returning the closure avoids redundant work per call.
- **Multi-start seed**: Fixed seed `np.random.default_rng(seed=0)` makes runs reproducible across restarts while still diversifying initialization.
- **Standard mode bounds approximation**: In standard (Re/Im) mode, bounding each Re and Im component to `[-a_max, a_max]` approximates `|w_n| ≤ a_max` but is not exact (the exact constraint is non-convex in Re/Im space). Documented in docstring; amplitude_only mode has exact bounds as an alternative.

**Open questions / known issues**:
- `STYLE.md` reshape example remains contradictory (not corrected to keep that file authoritative). The working code uses the empirically verified form.
- Bash tool failed on this platform (exit code 254, "stream closed before response"). All validation was run via PowerShell instead. Future scripts should prefer PowerShell on this machine.
- Phase wrapping: `np.angle()` returns the principal value in (−π, π]. A parsed phase of +191.825° stored as a complex phasor correctly round-trips; `np.angle()` returns −168.175° (differs by 2π). Not a bug, but could confuse manual inspection of the CSV output.
