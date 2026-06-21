# MATLAB Port

MATLAB port of the antenna array pattern optimization pipeline (`src/` +
`scripts/`). Built bottom-up and validated against the Python implementation via
the MATLAB MCP server — every core function has a unit test that compares its
output to a Python-generated reference fixture.

**Requires:** MATLAB R2026a + **Optimization Toolbox** (`fmincon`). No Signal
Processing Toolbox needed — the window functions are reimplemented from the
numpy/scipy definitions (`besseli` covers Kaiser).

## Entry points

```matlab
addpath('MATLAB/matlab_utils');          % from the repo root

run_optimization('config.yaml');         % full optimize pipeline -> results/optimizations/<ts>/
compare_classical('test_config.yaml');   % classical vs optimizer -> results/compare_classical/<ts>/
compare_classical('test_config.yaml', 8) % override n_side = 8
manual_weights_app('config.yaml');       % interactive weight-tuner GUI (blocks until window closed)
manual_weights_render('config.yaml', 'weights.csv', 'out.png');  % batch render from a weights CSV
manual_weights_render('config.yaml');    % uniform weights, default output path
```

The preferred way to run the scripts is via the entry-point wrappers in
`MATLAB/scripts/` (they call `addpath` for you):

```matlab
cd MATLAB/scripts
test_optimization_script   % drives run_optimization + compare_classical
manual_weights_app         % interactive GUI
```

Relative `element_patterns_dir` / `results_dir` in the config are resolved
against the repo root, so the scripts work regardless of the current folder.

## Layout (file = function, MATLAB convention)

All port functions live flat in `MATLAB/matlab_utils/`; `MATLAB/scripts/` holds
the thin entry-point wrappers (they `addpath('matlab_utils')` for you).

| Area | Files | Python origin |
|---|---|---|
| IO | `parse_cst_file`, `load_element_patterns` | `src/io/cst_parser.py` |
| Cost | `x_to_weights`, `weights_to_x`, `compute_array_factor`, `angular_window_mask`, `build_directive_physical_masks`, `build_cost_function`, `extended_grid_maps` | `src/cost/cost_function.py` |
| Optimize | `run_optimizer` (fmincon) | `src/optimize/optimizer.py` |
| Metrics | `evaluate_metrics`, `compute_directivity_dbi_grid`, `compute_hpbw`, `nearest_index` | `src/metrics/metrics.py` |
| Plot | `save_all_plots`, `save_pattern_gif`, `plot_comparison` | `src/plot/plotter.py` + compare_classical |
| compare_classical math | `build_ura_element_patterns`, `steering_phase_vector`, `data_driven_steering_vector`, `classical_weights`, `window_1d`, `principal_plane_cut`, `principal_plane_theta_axis`, `power_normalize_weights`, `evaluate_directive_metrics`, `run_scenario` | `scripts/compare_classical.py` |
| Shared | `read_config_yaml` (minimal YAML reader), `get_directive_field` | replaces PyYAML / `dict.get` |
| Scripts | `run_optimization`, `compare_classical`, `manual_weights_render` | `scripts/*.py` |
| Interactive GUI | `ManualWeightsTuner` (handle class), `scripts/manual_weights_app` | `scripts/manual_weights.py` |

## Tests

```matlab
runtests('MATLAB/tests')        % all 12 files, 33 tests
```

Each `tests/test_*.m` adds `MATLAB/matlab_utils` to the path, then loads a JSON fixture produced by
`tests/gen_reference_fixtures.py` (run with the Python project) and asserts the
MATLAB output matches within tolerance (exact for shapes/masks/indices,
~1e-9 for float math). Regenerate fixtures with:

```bash
python MATLAB/tests/gen_reference_fixtures.py            # all
python MATLAB/tests/gen_reference_fixtures.py optimizer  # one
```

## Conventions / differences from Python

- **Directives = cell array of structs** (`{struct(...), ...}`) — the analogue of
  Python's `list[dict]` with heterogeneous optional keys.
- **Element-pattern stack = `(N_elements, N_theta, N_phi)`**, matching the Python
  axis order; `compute_array_factor` sums over dim 1.
- **1-based indexing**: `nearest_index` / `run_optimizer.best_run_index` are
  1-based (Python returns 0-based); reshape/encoding adjusted accordingly.
- **`null_depth_db = []`** for peak directives (Python `None`); `metrics.json`
  writes these as `null`.
- **Optimizer is not bit-identical to scipy.** `fmincon` ≠ L-BFGS-B and MATLAB
  `rng(0)` ≠ numpy `default_rng(0)`, so iterates and random restarts differ.
  Both minimise the same `J(x)`; the converged cost matches scipy's global
  optimum on well-posed problems (verified to ~1e-4 in `test_optimizer`).
- **Interactive GUI** is provided by `ManualWeightsTuner` (`uifigure`-based handle
  class) launched via `scripts/manual_weights_app.m`. It is a full port of
  `manual_weights.py` — live heatmap, per-element sliders, directive table,
  polarisation selector, hover readout, and Load CSV/Config buttons. The
  non-interactive batch renderer `manual_weights_render` remains available for
  headless/scripted use.
