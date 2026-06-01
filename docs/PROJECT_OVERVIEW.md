# Antenna Array Pattern Optimization — Project Overview

## Purpose

This tool optimizes the complex element weights of a phased-array antenna to
shape its far-field radiation pattern according to user-defined beam and null
targets. It is intended for offline, pattern-in-the-loop design using measured
element patterns exported from CST Studio Suite.

The core algorithm is gradient descent (L-BFGS-B) on a composite cost function
built from per-directive power integrals. Weights are complex numbers — one per
array element — encoding both amplitude tapering and phase steering.

A future MATLAB port is planned; all Python-specific constructs are flagged with
`# [MATLAB]` inline comments.

---

## Scope

**In scope**

- Loading and parsing CST far-field exports (co-pol, cross-pol, total field).
- Formulating and minimizing a user-defined beam-shaping cost function.
- Multi-start L-BFGS-B optimization with amplitude and/or phase-only modes.
- Saving optimized weights (CSV), metrics (JSON), a run report, and plots.
- Interactive manual weight tuning with a live 2-D radiation pattern display.
- Benchmarking the optimizer against classical aperture tapering techniques
  (Hamming, Hanning, Kaiser, Chebyshev, Taylor) on a synthetic or real URA.

**Out of scope**

- Mutual coupling correction.
- Dynamic / real-time jammer tracking (planned as a separate project).
- Antenna synthesis or element placement optimization.

---

## Repository Layout

```
config.yaml               Main run configuration (edit this before each run)
scripts/
  run_optimization.py     Full pipeline entry point
  manual_weights.py       Interactive weight-tuning GUI
  compare_classical.py    Optimizer-vs-classical benchmark
  test_config.yaml        Config for compare_classical.py
src/
  io/cst_parser.py        CST .txt file parser
  cost/cost_function.py   Array factor + cost function
  optimize/optimizer.py   L-BFGS-B wrapper + multi-start
  plot/plotter.py         All matplotlib visualization
  metrics/metrics.py      Post-run per-directive scoring
data/element_patterns/    CST far-field exports (one .txt per element)
results/                  Timestamped output folders
docs/                     Design docs and session log
tests/                    pytest unit tests
```

Module boundaries are strict: no plotting in `optimize/`, no optimization
in `plot/`.

---

## Key Concepts

### Element patterns

The CST parser loads one ASCII `.txt` file per array element. Each file
contains an (N_theta × N_phi) grid of E-field values. The complex co-pol
pattern used in the optimizer is reconstructed as:

```
E_element(θ, φ) = Abs(Copol) · exp(j · Phase(Copol) · π/180)
```

All files in `element_patterns_dir` are loaded; the number of elements is
inferred from the file count.

### Array factor

The far-field of the full array is the coherent superposition of weighted
element patterns:

```
AF(θ, φ) = Σ_n  w_n · E_n(θ, φ)
```

where `w_n = a_n · exp(j·ψ_n)` is the complex weight for element n.

### Directives

Directives are the user-facing specification of what the optimizer should do.
Each directive defines a region of the (θ, φ) sphere and a goal:

| Field       | Description                                                    |
|-------------|----------------------------------------------------------------|
| `type`      | `"peak"` — maximize gain in the window; `"null"` — suppress it |
| `theta`     | Window center elevation (degrees, 0° = boresight)             |
| `phi`       | Window center azimuth (degrees)                               |
| `theta_width` / `phi_width` | Independent angular window extents (degrees) |
| `weight`    | Relative priority in the composite cost (higher = stronger)   |
| `aggregation` | How power is summarized: `"mean"` (default), `"max"`, `"min"` |

### Cost function

The optimizer minimizes a scalar cost `J`:

```
J(x) = Σ_k  λ_k · C_k(x)
```

Peak directive: `C_peak = −mean( |AF|² )` over the angular window (negative
because L-BFGS-B minimizes).

Null directive: `C_null = +mean( |AF|² )` over the window.

All internal computation is on **linear power**; dB conversion is display-only.

### Optimization

The solver is `scipy.optimize.minimize(method='L-BFGS-B')`. The variable
vector is the real/imaginary parts of all weights concatenated: length `2N`.

**Multi-start** is used to reduce the risk of local minima:

- Run 1 always starts from uniform weights (amplitude = 1, phase = 0).
- Runs 2…`n_restarts` start from random unit-amplitude phases.
- Optionally, one additional run per element starts with that element active
  and all others zeroed (`use_single_element_init`).
- All restarts use a fixed seed — results are fully reproducible.

The best result (lowest final cost J) across all restarts is kept.

### Polarization modes

| Mode    | Behavior                                                              |
|---------|-----------------------------------------------------------------------|
| `copol` | Optimize `AF = Σ w_n · E_copol_n`. Default.                          |
| `cross` | Optimize `AF = Σ w_n · E_cross_n`.                                   |
| `total` | Optimize incoherent sum `|AF_copol|² + |AF_cross|²`.                 |

### Physics conventions

- Internal computation: **radians**. Display and CSV output: **degrees**.
- Optimization operates on **linear field magnitude**; dB is for display only.
- Element 0 is the reference element (no per-element normalization).
- 0° elevation = boresight / zenith (z-axis).

---

## Outputs (per `run_optimization.py` run)

All outputs are saved to a timestamped subfolder under `results/`:

| File                  | Description                                        |
|-----------------------|----------------------------------------------------|
| `weights.csv`         | Optimized weights: index, amplitude, phase, Re, Im |
| `metrics.json`        | Per-directive and global gain metrics              |
| `run_report.txt`      | Human-readable run summary for cross-run comparison|
| `config.yaml`         | Snapshot of the config used for this run           |
| `pattern_polar.png`   | Polar pattern cut with directive overlays          |
| `pattern_cartesian.png` | Cartesian dB-vs-angle cut                        |
| `pattern_2d.png`      | Full (θ, φ) heatmap with directive markers         |
| `weights_amp.png`     | Amplitude bar chart per element                    |
| `weights_phase.png`   | Phase bar chart per element                        |
| `cost_history.png`    | Cost J vs iteration for all restarts               |
| `pattern_evolution.gif` | Optional: pattern animation over optimizer iters |

---

## Authoritative References

- [docs/pipeline.md](pipeline.md) — detailed data flow and stage-by-stage format spec
- [docs/STYLE.md](STYLE.md) — coding conventions and MATLAB porting flags
- [docs/notes.md](notes.md) — physics assumptions, open questions, session log
- [config.yaml](../config.yaml) — full annotated configuration schema
