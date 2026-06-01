# How to Use — Antenna Array Pattern Optimization

## Prerequisites

```
pip install numpy scipy matplotlib pyyaml pillow
```

Run all scripts from the **project root** (the folder containing `config.yaml`):

```
cd antenna_array_optimizations/
python scripts/<script_name>.py
```

---

## 1. `run_optimization.py` — Full Optimization Pipeline

Loads element patterns, runs L-BFGS-B optimization, and saves weights, metrics,
and plots to a timestamped folder under `results/`.

### Usage

```bash
python scripts/run_optimization.py                        # uses config.yaml
python scripts/run_optimization.py --config path/cfg.yaml # explicit config
```

### Quick configuration guide

Open `config.yaml` and set the four main sections before each run:

**Data paths**
```yaml
element_patterns_dir: "data/element_patterns/"  # folder with one CST .txt per element
```

**Directives** — what the optimizer should do

```yaml
directives:
  - type: "peak"
    theta: 30.0        # target elevation, degrees (0° = boresight)
    phi:   0.0         # target azimuth, degrees
    theta_width: 10.0  # angular window in elevation
    phi_width:  360.0  # angular window in azimuth (360° = rotationally symmetric)
    weight: 1.0        # relative cost priority

  - type: "null"
    theta: 60.0
    phi:   0.0
    theta_width: 5.0
    phi_width:  20.0
    weight: 5.0        # higher weight = stronger null enforcement
```

**Optimizer**
```yaml
optimizer:
  n_restarts:       5      # independent runs (more = better chance of global optimum)
  max_iterations:   500
  amplitude_bounds: [0.0, 1.0]
  phase_only:       false  # true = fix amplitudes to 1, optimize phases only
```

**Output** — toggle which plots to save; `results_dir` is the output root.

See `config.yaml` for the full annotated schema with all optional keys and defaults.

### Worked example: 30° steered beam with a null at 60°

```yaml
element_patterns_dir: "data/element_patterns/"
polarization: "copol"

directives:
  - type: "peak"
    theta: 30.0
    phi:   0.0
    theta_width: 15.0
    phi_width:  360.0
    weight: 1.0

  - type: "null"
    theta: 60.0
    phi:   0.0
    theta_width: 10.0
    phi_width:   30.0
    weight: 5.0

optimizer:
  n_restarts:       5
  max_iterations:   500
  cost_tolerance:   1.0e-8
  amplitude_bounds: [0.0, 1.0]
  phase_only:       false

output:
  results_dir: "results/"
  save_polar_plot:       true
  save_cartesian_plot:   true
  save_2d_projection_plot: true
  save_weight_plots:     true
  save_cost_history_plot: true
  save_pattern_gif:      false
```

Run and inspect the result:

```bash
python scripts/run_optimization.py
# Output: results/2026-05-27_143012/
#   weights.csv, metrics.json, run_report.txt, *.png
```

`run_report.txt` shows peak gain and null depth per directive in a human-readable
table — use this for quick cross-run comparison without opening plots.

---

## 2. `manual_weights.py` — Interactive Weight-Tuning GUI

Opens a tkinter window with a live 2-D (θ, φ) radiation pattern heatmap. Adjust
element weights interactively and see the pattern update in real time.

### Usage

```bash
python scripts/manual_weights.py                        # uses config.yaml
python scripts/manual_weights.py --config path/cfg.yaml
```

### GUI panels

| Panel | What it shows / lets you do |
|-------|-----------------------------|
| **2-D Radiation Pattern** | Heatmap of the array factor over the full sphere. Hover to read out gain at the cursor position. |
| **Element Weights** | Amplitude and phase entry + slider per element. Type a value or drag the slider — the pattern updates immediately. Click **Solo** to isolate one element. |
| **Directives** | Add/remove peak and null targets. Windows are overlaid on the heatmap (green = peak, red = null). Live gain/depth readout per directive. |
| **Metrics** | Global peak (dBi), peak angle, 3-dB HPBW, and total cost J. |

**Toolbar controls:**

- **Polarisation** — switch between `copol`, `cross`, and `total` (incoherent power sum).
- **Display** — `relative` (0 dB at live peak) or `absolute` (dBi with user-set color range).
- **Uniform Weights** — reset all weights to amplitude = 1, phase = 0.
- **Load Weights CSV** — load a `weights.csv` produced by `run_optimization.py` to
  visually inspect the optimized result, then fine-tune manually.

### Worked example: inspect optimized weights, then tweak the null manually

1. Run `run_optimization.py` to produce `results/2026-05-27_143012/weights.csv`.
2. Launch the GUI:
   ```bash
   python scripts/manual_weights.py
   ```
3. Click **Load Weights CSV** and open the file from step 1.
4. In the Directives panel, add a null at θ = 60°, φ = 0°, width = 10°.
   Observe the null depth in the inline readout.
5. Drag the phase slider of one element to see how the null depth responds.
6. The pattern heatmap and all metrics update on every change.

---

## 3. `compare_classical.py` — Optimizer vs. Classical Tapering Benchmark

Runs the L-BFGS-B optimizer alongside eight classical aperture tapering
techniques (Uniform, Hamming, Hanning, Kaiser β=3, Kaiser β=6, Chebyshev 25 dB,
Chebyshev 40 dB, Taylor 25 dB) on the same array and scenario(s), then plots a
side-by-side comparison.

Works with either a **synthetic** isotropic N×N URA (no real data needed) or
**real CST element patterns** loaded from a folder.

### Usage

```bash
python scripts/compare_classical.py                          # uses scripts/test_config.yaml
python scripts/compare_classical.py --config path/cfg.yaml  # explicit config
python scripts/compare_classical.py --n 8                   # override n_side (8×8 array)
```

The output figure is saved to `scripts/compare_classical_<N>x<N>.png` and displayed.

### Configuration (`scripts/test_config.yaml`)

Key sections:

**Array geometry**
```yaml
n_side: 4           # 4×4 = 16 elements (ignored when element_source: "folder")
d_over_lambda: 0.5  # half-wavelength spacing
```

**Element source**
```yaml
element_source: "folder"          # "synthetic" or "folder"
element_patterns_dir: "data/element_patterns/"
polarization: "copol"             # copol / cross / total
```

**Scenarios** — each scenario is one row in the comparison figure:
```yaml
scenarios:
  - label:           "Broadside beam, null at 20°"
    steer_theta_deg: 0.0
    steer_phi_deg:   0.0
    null_theta_deg:  20.0         # drawn as a vertical marker on the pattern plot
    directives:
      - type: "peak"
        theta: 0.0
        phi:   180.0
        theta_width: 5.0
        phi_width:   360.0
        weight: 1.0
      - type: "null"
        theta: 20.0
        phi:   0.0
        theta_width: 5.0
        phi_width:   360.0
        weight: 3.0
```

### Reading the output

The script prints a metrics table to the terminal for each scenario and directive:

```
  Technique          Peak (dBi)   Min* (dB)  Mean (dB)  Max (dB)
  ------------------  ----------  ----------  ----------  ----------
  uniform                  12.0       -3.2*       -1.8       0.0
  hamming                  10.5       -3.0*       -1.7       0.0
  ...
  optimized                11.8       -2.1*       -0.9       0.0 <--
```

For **peak** directives, `Min*` is the worst-case coverage point inside the window.
For **null** directives, `Max*` is the worst-case leakage inside the null region.

The figure shows:
- **Left column**: principal-plane pattern overlay (absolute dBi) for all techniques.
- **Right column**: whisker chart per directive — dot = solid-angle-weighted mean,
  whisker = [min, max], critical endpoint highlighted in orange/crimson.

### Worked example: compare a 4×4 synthetic URA for broadside peak + null at 60°

```yaml
# scripts/test_config.yaml
n_side: 4
d_over_lambda: 0.5
theta_step_deg: 1.0
phi_step_deg:   5.0
element_source: "synthetic"
n_restarts:     5
max_iterations: 500
cost_tolerance: 1.0e-9
plot_dynamic_range_db: 50

scenarios:
  - label:           "Broadside peak, null at 60°"
    steer_theta_deg: 0.0
    steer_phi_deg:   0.0
    null_theta_deg:  60.0
    directives:
      - type: "peak"
        theta: 0.0
        phi:   180.0
        theta_width: 10.0
        phi_width:   360.0
        weight: 1.0
      - type: "null"
        theta: 60.0
        phi:   0.0
        theta_width: 10.0
        phi_width:   30.0
        weight: 3.0
```

```bash
python scripts/compare_classical.py
# Prints metrics table to terminal
# Saves: scripts/compare_classical_4x4.png
```

---

## Typical Workflow

```
1. Export element far-field patterns from CST → data/element_patterns/

2. Edit config.yaml:
   - Set element_patterns_dir
   - Define peak / null directives
   - Set optimizer.n_restarts (5–10 for production)

3. python scripts/run_optimization.py
   → results/<timestamp>/weights.csv  (optimized weights)
   → results/<timestamp>/run_report.txt  (quick gain / null-depth summary)

4. python scripts/manual_weights.py
   → Load the weights.csv from step 3
   → Inspect the live 2-D pattern, check directive overlays

5. (Optional) python scripts/compare_classical.py
   → Quantify how much the optimizer improves on classical tapering
```
