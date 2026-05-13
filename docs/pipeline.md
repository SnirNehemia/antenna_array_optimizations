# Antenna Array Pattern Optimization — Pipeline

## Overview

This document describes the full data and processing pipeline, from raw CST exports to
optimized complex weights and evaluated radiation patterns.

---

## Directory Structure

```
antenna_array_optimizer/
│
├── CLAUDE.md           # Project context for Claude Code (root level)
├── config.yaml         # User-facing run configuration
│
├── src/
│   ├── io/             # CST file parsing and data ingestion
│   │   └── cst_parser.py
│   ├── cost/           # Objective function and target mask construction
│   │   └── cost_function.py
│   ├── optimize/       # Gradient descent loop and weight solver
│   │   └── optimizer.py
│   ├── plot/           # All visualization routines
│   │   └── plotter.py
│   └── metrics/        # Post-run evaluation and per-directive scoring
│       └── metrics.py
│
├── data/
│   └── element_patterns/   # Raw CST exports (.txt / .ffs), one file per element
│
├── tests/              # Unit tests (pytest)
├── docs/               # This file, notes.md, STYLE.md, and any design docs
├── scripts/            # Top-level runner scripts / entry points
│   └── run_optimization.py
└── results/            # Timestamped output folders (weights, plots, metrics)
```

---

## Stage 0 — CST Export Format

CST Studio exports element far-field patterns as ASCII `.txt` files.

### Filename convention (CST default)
```
farfield__f_<freq>__<element_index>__E_field_<res>degsres.txt
```

### File structure
- **Row 1**: Column headers (tab/space-separated)
- **Row 2**: Separator line (`---...`)
- **Rows 3+**: Data, one sample per line, 8 columns:

| Column | Name            | Unit   | Description                          |
|--------|-----------------|--------|--------------------------------------|
| 0      | Theta           | deg    | Elevation angle (0° = boresight)     |
| 1      | Phi             | deg    | Azimuth angle                        |
| 2      | Abs(E)          | V/m    | Total E-field magnitude              |
| 3      | Abs(Cross)      | V/m    | Cross-polarization magnitude         |
| 4      | Phase(Cross)    | deg    | Cross-polarization phase             |
| 5      | Abs(Copol)      | V/m    | Co-polarization magnitude            |
| 6      | Phase(Copol)    | deg    | Co-polarization phase                |
| 7      | Ax.Ratio        | —      | Axial ratio                          |

### Angular grid (example: 5° resolution)
- Theta: 0° to 180°, step = 5° → 37 unique values
- Phi:   0° to 355°, step = 5° → 72 unique values
- Total samples: 37 × 72 = 2,664 rows

> The parser must be resolution-agnostic: it auto-detects the step size from
> the unique sorted values of Theta and Phi, so it works for any export resolution.

### Which field to use for optimization
By default, use **Abs(Copol)** and **Phase(Copol)** to construct the complex
element pattern used in the array factor computation:

```
E_element(theta, phi) = Abs(Copol) * exp(j * Phase(Copol) * pi/180)
```

This can be switched to total-field `Abs(E)` via a config flag. Cross-pol
optimization is a future extension.

---

## Stage 1 — IO / Parsing (`src/io/cst_parser.py`)

**Input**: Path to one or more CST `.txt` files (one per array element).

**Process**:
1. Skip the 2-line header.
2. Parse all 8 columns as float.
3. Detect angular resolution from unique Theta and Phi values.
4. Reshape flat data into a 2D grid: shape `(N_theta, N_phi)` for each column.
5. Construct the complex element pattern grid from Copol magnitude + phase.

**Output** (per element):
```python
{
    "theta_deg":   np.ndarray, shape (N_theta,)          # elevation angles
    "phi_deg":     np.ndarray, shape (N_phi,)             # azimuth angles
    "E_abs":       np.ndarray, shape (N_theta, N_phi)     # total magnitude [V/m]
    "copol_abs":   np.ndarray, shape (N_theta, N_phi)     # copol magnitude [V/m]
    "copol_phase": np.ndarray, shape (N_theta, N_phi)     # copol phase [deg]
    "E_complex":   np.ndarray, shape (N_theta, N_phi)     # copol complex field
}
```

---

## Stage 2 — Array Model

For an N-element array, the total array pattern is the coherent sum of weighted
element patterns:

```
AF(theta, phi) = sum_{n=0}^{N-1}  w_n * E_n(theta, phi)
```

where:
- `w_n = a_n * exp(j * psi_n)` — complex weight for element n
- `a_n` — amplitude (real, ≥ 0)
- `psi_n` — phase shift [radians]
- `E_n(theta, phi)` — complex element pattern from Stage 1

The optimization variable vector is the concatenated real/imaginary parts of all
weights: `x = [Re(w_0), Im(w_0), Re(w_1), Im(w_1), ..., Re(w_{N-1}), Im(w_{N-1})]`
(length `2N`), which is compatible with gradient-based solvers like L-BFGS-B.

---

## Stage 3 — Cost Function (`src/cost/cost_function.py`)

The optimizer minimizes a scalar cost `J(x)`:

```
J(x) = sum_k  lambda_k * C_k(x)
```

where each directive k contributes a term `C_k`:

### Peak directive
Maximize gain at target angle `(theta_t, phi_t)` over width `delta`:
```
C_peak = -mean( |AF(theta, phi)|^2  for (theta,phi) in angular window )
```
(negative because L-BFGS-B minimizes)

### Null directive
Suppress gain at target angle `(theta_t, phi_t)` over width `delta`:
```
C_null = mean( |AF(theta, phi)|^2  for (theta,phi) in angular window )
```

### Composite cost
```
J(x) = sum_k  lambda_k * C_k(x)
```
where `lambda_k` is the user-supplied `weight` for directive k (default = 1.0).

All computations are in **linear field magnitude**; dB conversion is only for display.

---

## Stage 4 — Optimizer (`src/optimize/optimizer.py`)

**Solver**: `scipy.optimize.minimize` with `method='L-BFGS-B'`

**Variables**: `x` — real vector of length `2N` (Re/Im of each weight)

**Optional constraints**:
- Amplitude bounds per element (e.g., `0 ≤ a_n ≤ 1`)
- Phase-only mode (fix amplitudes to 1, optimize phases only)

**Initialization**: By default, uniform weights (all ones, zero phase). Optionally
random multi-start to escape local minima.

**Output**:
- `weights_complex`: `np.ndarray`, shape `(N,)`, complex — final optimized weights
- `cost_history`: list of cost values per iteration
- `result`: full `scipy.OptimizeResult` object

---

## Stage 5 — Metrics (`src/metrics/metrics.py`)

After optimization, evaluate per-directive performance:

| Directive | Metric                                    | Unit |
|-----------|-------------------------------------------|------|
| peak      | Realized gain at target angle             | dBi  |
| null      | Null depth at target angle                | dB   |
| global    | Peak-to-null ratio (if both are defined)  | dB   |

Also report: total cost `J`, per-directive cost breakdown, and iteration count.

---

## Stage 6 — Visualization (`src/plot/plotter.py`)

All plots are generated by functions in `plotter.py`. Each function saves to a
timestamped folder under `results/`.

| Plot                    | Format                  | Description                              |
|-------------------------|-------------------------|------------------------------------------|
| Array pattern (polar)   | Polar, dB scale         | |AF|² vs angle, overlaid with target mask |
| Array pattern (Cartesian)| dB vs angle (cut)      | Theta-cut at fixed Phi (default Phi=0°)  |
| Weights — amplitude     | Bar chart               | |w_n| per element                         |
| Weights — phase         | Bar chart               | ∠w_n in degrees per element              |
| Cost history            | Line plot               | J vs iteration                           |

---

## Stage 7 — Entry Point (`scripts/run_optimization.py`)

The top-level runner:
1. Loads config (element pattern paths, optimization directives, solver settings).
2. Calls `cst_parser` → loads element patterns.
3. Calls `optimizer` → runs gradient descent.
4. Calls `metrics` → evaluates result.
5. Calls `plotter` → saves all plots.
6. Writes `results/<timestamp>/weights.csv` and `results/<timestamp>/metrics.json`.

---

## Data Flow Diagram

```
CST .txt files
     │
     ▼
[Stage 1: cst_parser]
     │  E_complex[N_theta, N_phi] per element
     ▼
[Stage 2: Array Model]
     │  AF(theta, phi; x)
     ▼
[Stage 3: cost_function] ◄── directives (peaks, nulls, weights)
     │  J(x)
     ▼
[Stage 4: optimizer] ──────── scipy L-BFGS-B
     │  weights_complex[N]
     ▼
[Stage 5: metrics] ──────────────────────────────► metrics.json
     │
[Stage 6: plotter] ──────────────────────────────► plots/
```
