# Antenna Array Pattern Optimization

Gradient-descent beam shaping for phased arrays from CST Studio far-field exports.
Define peak and null directives in a YAML config; the optimizer finds the complex
element weights that satisfy them using L-BFGS-B with multi-start.

---

## What it does

- Parses CST far-field exports (co-pol, cross-pol, or total field) for each array element.
- Minimizes a composite cost function built from user-defined peak and null windows.
- Saves optimized weights (CSV), gain/null metrics (JSON), a run report, and plots.
- Provides an interactive GUI to inspect or manually tune weights with a live 2-D pattern.
- Benchmarks the optimizer against classical aperture tapering techniques (Hamming,
  Kaiser, Chebyshev, Taylor) on synthetic or real element patterns.

---

## Setup

```bash
pip install numpy scipy matplotlib pyyaml pillow
```

Run all scripts from the project root (the folder containing `config.yaml`).

---

## Quick start

**1. Run the optimizer**
```bash
# Edit config.yaml first: set element_patterns_dir and define directives
python scripts/run_optimization.py
# Results saved to results/<timestamp>/
```

**2. Inspect results interactively**
```bash
python scripts/manual_weights.py
# Load the weights.csv from the results folder via the toolbar button
```

**3. Compare against classical tapering**
```bash
# Edit scripts/test_config.yaml to set array size, scenarios, and element source
python scripts/compare_classical.py
```

---

## Documentation

| Document | Content |
|----------|---------|
| [docs/PROJECT_OVERVIEW.md](docs/PROJECT_OVERVIEW.md) | Scope, key concepts (array factor, directives, cost function), architecture, and output file reference |
| [docs/HOW_TO_USE.md](docs/HOW_TO_USE.md) | Step-by-step guide for all three scripts, config key reference, and worked examples |
| [docs/OFFLINE_SETUP.md](docs/OFFLINE_SETUP.md) | Installing on an air-gapped / offline machine (wheel download + transfer instructions) |
| [docs/pipeline.md](docs/pipeline.md) | Detailed data flow and per-stage format specification |
| [config.yaml](config.yaml) | Fully annotated run configuration schema |

---

## MATLAB port

A full MATLAB port lives in [`MATLAB/`](MATLAB/README.md). It replicates the
entire pipeline — CST parser, cost function, optimizer, metrics, plots,
compare-classical benchmark, and the interactive weight-tuner GUI — using
only built-in MATLAB functions plus the Optimization Toolbox.

```matlab
addpath('MATLAB');
run_optimization('config.yaml');          % optimize → results/optimizations/<ts>/
compare_classical('test_config.yaml');    % benchmark vs classical tapering
manual_weights_app('config.yaml');        % interactive weight-tuner GUI
```

See [`MATLAB/README.md`](MATLAB/README.md) for requirements, layout, and test instructions.

---

## Repository layout

```
config.yaml               Main run configuration
scripts/
  run_optimization.py     Full pipeline entry point
  manual_weights.py       Interactive weight-tuning GUI
  compare_classical.py    Optimizer-vs-classical benchmark
  test_config.yaml        Config for compare_classical.py
src/
  io/                     CST parser
  cost/                   Cost function
  optimize/               L-BFGS-B wrapper + multi-start
  plot/                   Visualization
  metrics/                Post-run scoring
data/element_patterns/    CST far-field exports (one .txt per element)
results/                  Timestamped output folders
MATLAB/                   Complete MATLAB port (see MATLAB/README.md)
```
