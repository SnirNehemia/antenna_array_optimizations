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

### 2026-06-01 — Fix evaluate_metrics for pole-crossing directives

**Problem**: For a null directive at θ=−30°, φ=0° (back hemisphere via pole),
the run report and GUI inline label showed `gain_dbi = +11.47 dBi` and
`null_depth = −4.80 dB`, while the displayed red-box max was ~−8 dBi. The null
was actually achieved (~24 dB deep); only the metric was wrong.

**Root cause** (`src/metrics/metrics.py`, `evaluate_metrics`): the metric sampled
the directive on the **raw** physical θ∈[0,180°] grid, unlike the cost function and
the GUI overlay which both use the extended grid. `_nearest_index(θ, −30)` clamped
to θ=0° (boresight, next to the main beam) → bogus +11 dBi. The window mask
(`angular_window_mask` on the raw grid, window θ∈[−35,−25]) was entirely off-grid →
empty → `mean_window_power = 0` → `cost_term = 0.0` and a meaningless null depth.

**Fix**: `evaluate_metrics` now builds masks via `build_directive_physical_masks`
(the same extended-grid masks the optimizer and overlay use), so θ=−30°,φ=0° maps
to its physical mirror θ=30°,φ=180°. `gain_dbi` is reported as the **max
directivity inside the physical window** (peak → achieved gain; null → worst-case
leakage), matching the brightest pixel of the on-screen box. Window mean power for
`cost_term` uses the same mask. Verified against the
`2026-06-01_094337` run: null now reads −7.90 dBi / depth −23.95 dB / cost_term 19.1,
peak 16.06 dBi, peak-to-null 23.95 dB.

**Decisions made**:
- `gain_dbi` semantics changed from "value at nearest center grid point" to
  "max inside the angular window" for **both** peak and null, so every inline
  directive label matches the brightest pixel of its on-screen window box.
  Kept a nearest-mapped-center fallback only for a degenerate (empty) window.
- `cost_term` still uses solid-angle-weighted mean (the default `"mean"`
  aggregation); per-directive `aggregation: max/min` is still not reflected in the
  reported `cost_term` (pre-existing limitation, unchanged).

### 2026-05-27 — Extended-grid angle wrap-around, aggregation mode, theta-axis flip, physical mask overlays

**Implemented**:

- **Phi 0°/360° and theta pole-crossing wrap-around** (`src/cost/cost_function.py`):
  Extended-grid approach mirrors the physical pattern at both poles (θ<0°, θ>180°)
  and tiles phi three times. Masks are computed on the extended grid with a plain
  `abs(angle − target) ≤ hw` comparison — no special-case code needed.
  Physical-index maps (`theta_phys_idx`, `phi_offset`) are precomputed once;
  `cost_fn` samples `power_grid` at sparse mask-True points only
  (`power_grid[phys_theta, phys_phi]`), so per-iteration overhead is O(K) where K
  is the number of in-window samples rather than O(9·N_theta·N_phi).
  Phi offset for mirrored theta rows = `N_phi // 2` (180° shift, valid for the
  standard 1°-step 360-point phi grid).
  Theta pole-crossing known limitation: wrap is handled correctly but requires both
  theta halves to be distinct grid entries; no special handling for directives
  exactly at the poles.

- **`aggregation` per-directive config key** (`src/cost/cost_function.py`, `config.yaml`):
  `"mean"` (default, solid-angle-weighted, current behaviour), `"max"` (worst-case
  point — best for null suppression), `"min"` (best-case point — best for flat-beam
  enforcement). Validated at `build_cost_function` call time.

- **`build_directive_physical_masks`** (`src/cost/cost_function.py`):
  Standalone public function returning one `(N_theta, N_phi)` bool array per
  directive, built with the same extended-grid logic as `build_cost_function`.
  Used by the plotter and manual tuner to display exactly the area the optimizer
  sees, including wrapped regions.

- **2D heatmap theta axis flipped** (`src/plot/plotter.py`): `y_lim = (180.0, 0.0)`
  so θ=0° (boresight) is at the top in both the static PNG and the GIF animation.
  Equal-area mode was already correct.

- **Physical mask overlays replace rectangles** (`src/plot/plotter.py`,
  `scripts/manual_weights.py`): `save_2d_projection_plot` and `save_pattern_gif`
  now call `build_directive_physical_masks` and draw `contourf` fill +
  `contour` border on the physical grid. The manual tuner (`_update_directive_overlays`)
  uses `pcolormesh` for the fill (reliable `.remove()` across all matplotlib versions)
  and `contour` for the border with a `.collections` fallback for pre-3.8 matplotlib.
  All three views show the same wrapped/pole-crossing window that the optimizer uses.

- **Run report `width` fix** (`scripts/run_optimization.py`): report line now reads
  `theta_width` / `phi_width` from the directive instead of the removed `width` key,
  so configs using per-axis widths no longer crash.

**Decisions made**:
- Extended grid is always built (not conditionally when wrapping is needed) — keeps
  code uniform and the precomputation cost is negligible (~microseconds).
- `aggregation: "max"` does not use solid-angle weighting (raw max over mask points)
  because the intent is to drive down the worst-case sidelobe regardless of its
  angular area. Same rationale for `"min"`.
- Theta pole-crossing wraps the physical phi to `phi + 180°` using integer index
  shift (`N_phi // 2`), which is exact for a uniform 1° grid; a comment marks this
  assumption for future datasets with different resolutions.

### 2026-05-26 — compare_classical.py: benchmark script, CST steering fix, dBi display

**Implemented**:
- `scripts/compare_classical.py` — new standalone benchmark script. Compares eight
  classical tapering/steering techniques (uniform, Hamming, Hanning, Kaiser β=3/6,
  Chebyshev 25/40 dB, Taylor 25 dB) against the L-BFGS-B optimizer on a configurable
  N×N URA. Outputs a two-column figure per scenario: left = principal-plane pattern
  overlay in absolute dBi, right = per-directive whisker chart showing min/mean/max
  gain inside each directive window.
- `scripts/test_config.yaml` — companion config file. Every test variable is exposed:
  `n_side`, `d_over_lambda`, angular grid steps, `element_source` (`"synthetic"` or
  `"folder"`), `element_patterns_dir`, `polarization`, optimizer settings
  (`n_restarts`, `max_iterations`, `cost_tolerance`), `plot_dynamic_range_db`, and a
  `scenarios` list (each scenario specifies `steer_theta_deg`, `steer_phi_deg`,
  `null_theta_deg`, and a `directives` list with the same schema as `config.yaml`).
- **Two element-source modes**:
  - `"synthetic"` — builds ideal isotropic URA phase-factor patterns from `n_side`
    and `d_over_lambda`. No real data required.
  - `"folder"` — loads CST Studio far-field exports via `src/io/cst_parser.py`.
    Infers `n_side` from `sqrt(n_elements)`. Theta/phi grids come from the files;
    `theta_step_deg`/`phi_step_deg` are ignored. Field component selected by
    `polarization` key.
- **Element pattern normalization** (`_load_folder_element_patterns`): after stacking,
  all patterns are divided by the global peak amplitude. CST exports have
  simulation-dependent absolute V/m amplitudes; without this step a 4×4 array would
  show ~30 dBi instead of the correct ~12 dBi bound for 16 isotropic-equivalent
  elements. The normalization preserves all relative phase and inter-element amplitude
  information.
- **Data-driven steering** (`_data_driven_steering_vector`): evaluates each element's
  pattern at the target direction and conjugates the phase. Used instead of the
  geometric `_steering_phase_vector` whenever real element patterns are provided. CST
  exports are phased relative to each element's own feed, not the array centre, so the
  geometric formula gives no inter-element progressive phase — classical techniques
  produced an unsteered broadside beam until this fix. Data-driven steering works for
  both synthetic and real patterns.
- **Absolute-dBi overlay plot** (`_plot_pattern_overlay`): all techniques plotted on a
  shared dBi y-axis so gain loss from tapering is immediately visible (e.g., Hamming
  at 7.6 dBi vs uniform at 12 dBi). Y-range: `[max_peak − dynamic_range_db,
  max_peak + 2]`. A −3 dB reference line is drawn.
- **Right-panel whisker charts** (`_plot_directive_whiskers`): y-axis changed from
  "dB relative to pattern peak" to absolute "Gain (dBi)" by adding each technique's
  `peak_dbi` back to the stored `min_db`/`mean_db`/`max_db` values. A dashed
  reference line marks the highest `peak_dbi` across all techniques. Critical-tip
  markers (▼ for peak-directive min, ▲ for null-directive max) now carry a small
  numerical label (1 decimal place in dBi) placed just below/above the marker.

**Decisions made**:
- Element pattern normalization is applied only in the folder-loading path, not to
  synthetic patterns (which are already unit-amplitude). This keeps the synthetic path
  as a clean mathematical reference.
- Data-driven steering is used unconditionally whenever element patterns are provided
  (both synthetic and folder). For synthetic patterns it yields identical results to
  the geometric formula, so there is no regression.
- The benchmark optimizer config always sets `use_single_element_init: False` to keep
  run time predictable; users control `n_restarts` from `test_config.yaml`.
- Critical-marker text is colored the same as the marker (darkorange for classical,
  crimson for optimized) for visual grouping.

**Open questions / known issues**:
- The "dBi" reported after element-pattern normalization is referenced to the
  simulation's peak-amplitude element, not a true isotropic radiator with 1 W input.
  For relative technique comparison this is self-consistent; absolute gain claims
  against a calibrated reference would need the CST patterns to carry a known input
  power normalization.

### 2026-05-26 — manual_weights.py GUI enhancements (hover, HPBW, axis invert)

**Implemented**:
- `src/metrics/metrics.py`:
  - `evaluate_metrics` now returns four new keys: `global_peak_theta_deg`,
    `global_peak_phi_deg` (location of the global peak), `hpbw_theta_deg`,
    `hpbw_phi_deg` (3 dB half-power beamwidth in the θ-cut and φ-cut at the peak).
  - `_compute_hpbw` extended with an optional `opposite_cut` parameter. When the
    left (or right) scan reaches the grid boundary without finding a 3 dB crossing,
    the search continues in the θ-cut at φ+180° (the other side of the pole). The
    virtual crossing index is set to `−opp_extra` (north pole) or `(n−1)+opp_extra`
    (south pole), so `right − left` gives the correct full beamwidth. Fixes HPBW for
    beams pointing near θ=0° or θ=180°. φ-cut wrap-around handled via `np.roll` (unchanged).
- `scripts/manual_weights.py`:
  1. **In-axes cursor annotation**: a semi-transparent text box overlaid in the top-left
     corner of the heatmap shows `θ=X° φ=Y° D=Z dB(i)` as the mouse moves. Value is in
     the active display mode's units (dB relative or dBi absolute). Hidden when the mouse
     leaves the axes (`axes_leave_event`).
  2. **θ-axis inversion**: `set_ylim(180.0, 0.0)` so 0° (boresight/zenith) is at the top
     (standard antenna convention).
  3. **Power-normalized weights**: before computing the array factor, weights are divided by
     `||w||₂` to match the optimizer's `cost_fn` convention; Total J in the GUI now equals
     the optimizer's objective. Directivity is scale-invariant so the displayed pattern is
     unchanged.
  4. **Metrics panel**: removed "Peak-to-null" row; added "Peak angle" (θ,φ of global peak)
     and "3 dB HPBW" (θ/φ beamwidths from `evaluate_metrics`).
  5. **Directive inline results**: each directive row now shows its live gain (for peaks) or
     `gain (null_depth_db)` (for nulls) in green/red next to the × button.
  6. **Status bar**: relocated inside the "2-D Radiation Pattern" LabelFrame (below canvas).

**Decisions made**:
- Power-normalization in the GUI: directivity `D = 4π|AF|²/P_total` is invariant to any
  overall amplitude scaling, so normalizing weights does not alter any visual output. The
  change only affects the Total J value, making it numerically consistent with the optimizer.
- Hover annotation uses the display-mode grid (`_last_display_grid`) so the shown value
  matches the colorbar exactly. The absolute dBi grid (`_dbi_grid`) is computed in parallel
  for metrics but not shown separately in the hover.
- `_compute_hpbw` returns 0 when the beam never drops 3 dB within the grid (e.g., isotropic
  radiator or omnidirectional pattern); this is correct — callers treat it as "HPBW > grid extent".

**Open questions / known issues**:
- None.

### 2026-05-19 — Solid-angle weighting in cost/metrics; flat-θ visualization

**Implemented**:
- `src/cost/cost_function.py` — `_directive_cost` now uses a solid-angle-weighted mean
  `Σ(|AF|² × sin θ) / Σ(sin θ)` over the directive mask instead of a uniform pixel mean.
  `sin_theta` is pre-computed once per `build_cost_function` call (outside the closure).
  `power_grid = |AF|²` computed once outside the directive loop per `cost_fn` call.
  A 5° window at θ=2° (pole) has the same pixel count as at θ=90° (equator) but ~29× less
  solid angle — the new weighting correctly de-emphasizes polar pixels proportionally.
- `src/metrics/metrics.py` — `evaluate_metrics` window mean uses the same sin θ weighting
  so `cost_term` reconstruction stays consistent with the optimizer.
- `src/plot/plotter.py` — `save_2d_projection_plot` default changed to `equal_area=False`
  (flat linear θ axis, 0–180°). The cos(θ) path is retained and still selectable via
  `equal_area=True` or `plot_equal_area: true` in config.
- `src/plot/plotter.py` — `save_pattern_gif` reverted to flat θ y-axis and flat-θ directive
  rectangle coordinates.
- `scripts/manual_weights.py` — `_build_pattern_panel` uses `self._theta_deg` as y-coordinates
  (reverted from cos θ); directive overlays in `_update_directive_overlays` draw in θ space.
  `self._cos_theta` attribute removed.
- `config.yaml` — `plot_equal_area: false`.

**Decisions made**:
- Flat θ axis is preferred for visualization: the cos(θ) equal-area projection compresses
  directive windows near the poles to near-zero visual height, which is confusing for
  pole-directed antenna optimization. Equal-area property is now captured correctly in the
  cost function (sin θ weighting) rather than in the visual projection.
- `gain_dbi` (directivity point lookup at target grid index) is unaffected — it is correct
  at all angles regardless of weighting.
- A directive exactly at θ=0° has zero solid-angle weight (sin 0° = 0); the optimizer will
  treat it as a zero-cost term. Physically correct (the pole is a single degenerate direction)
  but may surprise users — note in config if pole directives are needed.

**Open questions / known issues**:
- None new.

### 2026-05-18 — Cross-pol, absolute-dBi, CST directivity convention, separate θ/φ widths, GIF

**Implemented**:
- `src/io/cst_parser.py` — `parse_cst_file` extended to extract cross-polarization columns
  and return `cross_complex = cross_abs × exp(j × cross_phase_rad)`.
- `scripts/manual_weights.py` — polarization toolbar combobox extended with `cross` and
  `total` modes. `total` computes orthogonal power sum `|AF_copol|² + |AF_xpol|²`.
- `scripts/manual_weights.py` — display-mode toolbar added: `relative` (peak-normalized) and
  `absolute` (dBi, user-set clim via min/max entry widgets).
- `src/metrics/metrics.py` — `_compute_directivity_dbi_grid` gains optional `normalizer_power`
  kwarg; when provided it replaces the single-pol integral as denominator (CST/IEEE Std 149
  partial directivity convention). `evaluate_metrics` accepts `normalizer_power` and
  `precomputed_array_factor`.
- `scripts/manual_weights.py` — always computes both copol and cross AFs; passes
  P_copol + P_cross as `normalizer_power` to metrics and directivity. This ensures
  D_total ≥ D_copol ≥ 0 at all angles (fixed previous incorrect total < copol result).
- `src/cost/cost_function.py` — `angular_window_mask` signature changed from `width_deg` to
  `theta_width_deg, phi_width_deg`; `_directive_cost` call site updated. Directives support
  independent elevation and azimuth window widths; `width` remains as symmetric shorthand.
- `src/metrics/metrics.py`, `src/plot/plotter.py` — updated all `angular_window_mask` call
  sites to use `theta_width_deg`/`phi_width_deg` with `width` fallback.
- `scripts/manual_weights.py` — directive table now has `θW(°)` and `φW(°)` columns (split
  from single `W(°)` column).
- `src/optimize/optimizer.py` — `_run_single_optimization` callback records `xk.copy()` to
  `xk_history`; `run_optimizer` returns `weights_history` (decoded complex weights per best-run
  iteration) and `all_run_labels`, `all_cost_histories`, `best_run_index`.
- `src/plot/plotter.py` — new `save_pattern_gif()`: animated GIF showing radiation pattern
  evolving over optimizer iterations. Two-panel layout: 2D heatmap with directive overlays
  (left) + convergence plot with current-iteration dot (right). Auto-strides to
  `gif_max_frames`. Requires `pillow`.
- `scripts/run_optimization.py` — calls `save_pattern_gif` conditionally when
  `save_pattern_gif: true` in config output section. Import updated accordingly.
- `config.yaml` — added `theta_width`/`phi_width` fields to directive schema comments;
  added `plot_equal_area`, `save_pattern_gif`, `gif_max_frames` to output section.

**Decisions made**:
- CST partial directivity convention (IEEE Std 149): D_copol = 4π|AF_copol|² / (P_copol + P_cross).
  Ensures D_total ≥ D_copol always. `run_optimization.py` (single stack) uses the single-pol
  denominator as a documented approximation.
- `cross` polarization selectable in `run_optimization.py` via `polarization: "cross"` config;
  `total` is GUI-only (`run_optimization.py` only supports coherent single-stack optimization).
- GIF always uses flat θ axis (after the equal-area projection was initially used then reverted
  per the 2026-05-19 entry).

**Open questions / known issues**:
- GIF file size can be large for many iterations or high-resolution grids; `gif_max_frames`
  caps frame count but not file size.

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
