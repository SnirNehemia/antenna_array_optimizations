# ══════════════════════════════════════════════════════════════════
# RUN_OPTIMIZATION
# Top-level entry point: load config, run the full optimization pipeline,
# save weights, metrics, and plots to a timestamped results folder.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

# ────────────────────────── IMPORTS ───────────────────────────────
import argparse
import csv
import json
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

# Ensure the project root (parent of scripts/) is on sys.path so that
# `src.*` imports resolve when the script is invoked as a file rather
# than as a module (i.e. `python scripts/run_optimization.py`).
sys.path.insert(0, str(Path(__file__).parent.parent))
# [MATLAB] addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

import numpy as np
import yaml

from src.io.cst_parser import load_element_patterns
from src.cost.cost_function import compute_array_factor
from src.optimize.optimizer import run_optimizer
from src.metrics.metrics import evaluate_metrics
from src.plot.plotter import save_all_plots, save_pattern_gif

# ────────────────────────── CONSTANTS ─────────────────────────────

# Polarization modes wired through the optimizer pipeline.
# "copol" → E_complex (co-pol), "cross" → cross_complex (cross-pol).
# "total" is supported only in the interactive manual_weights tool, where the
# power-sum |AF_copol|² + |AF_xpol|² is rendered without a single coherent AF.
SUPPORTED_POLARIZATIONS = ("copol", "cross")

# Mapping from polarization name → pattern dict key used to build the stack.
POLARIZATION_TO_PATTERN_KEY = {
    "copol": "E_complex",
    "cross": "cross_complex",
}

# Required top-level keys in config.yaml.
REQUIRED_CONFIG_KEYS = ("element_patterns_dir", "directives", "optimizer", "output")


# ────────────────────────── HELPERS ───────────────────────────────

def _load_config(config_path):
    """Load config.yaml and return its contents as a dict.

    Args:
        config_path (str | Path): Path to the YAML configuration file.

    Returns:
        dict: Parsed configuration dictionary.

    Raises:
        FileNotFoundError: If the config file does not exist at config_path.
    """
    with open(config_path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)
    # [MATLAB] config = yaml_load(config_path);  % requires a YAML toolbox


def _validate_config(config):
    """Check that all required config keys are present and polarization is supported.

    Args:
        config (dict): Parsed configuration dictionary.

    Raises:
        KeyError: If a required top-level key is absent.
        ValueError: If the specified polarization is not yet implemented.
    """
    for key in REQUIRED_CONFIG_KEYS:
        if key not in config:
            raise KeyError(
                f"Missing required config key: '{key}'. "
                f"Required keys: {REQUIRED_CONFIG_KEYS}."
            )

    polarization = config.get("polarization", "copol")
    if polarization not in SUPPORTED_POLARIZATIONS:
        raise ValueError(
            f"Polarization '{polarization}' is not yet implemented. "
            f"Supported options: {SUPPORTED_POLARIZATIONS}. "
            "Set polarization: \"copol\" in config.yaml."
        )


def _create_output_dir(results_dir):
    """Create and return a timestamped output subfolder under results_dir.

    Args:
        results_dir (str | Path): Base results directory from config.

    Returns:
        Path: Newly created timestamped output directory.
    """
    timestamp  = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    # [MATLAB] timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
    output_dir = Path(results_dir) / timestamp
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _save_weights_csv(weights_complex, output_path):
    """Save optimized complex weights to a CSV file.

    Columns: element_index, amplitude, phase_deg, real, imag.

    Args:
        weights_complex (np.ndarray): Optimized complex weights, shape (N_elements,).
            Units: dimensionless (V/V).
        output_path (str | Path): Full path for the CSV output file.
    """
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["element_index", "amplitude", "phase_deg", "real", "imag"])
        for n, w in enumerate(weights_complex):
            writer.writerow([
            n,
            f"{np.abs(w):.2f}",
            f"{np.angle(w, deg=True):.2f}",
            f"{w.real:.2f}",
            f"{w.imag:.2f}"
        ])
    # [MATLAB] writematrix([indices, amps, phases, reals, imags], output_path);


def _save_metrics_json(metrics, output_path):
    """Save the metrics dict to a JSON file with 2-space indentation.

    All values in the metrics dict are native Python types (float, int, str, None),
    so no custom JSON encoder is required.

    Args:
        metrics (dict): Output of evaluate_metrics(...).
        output_path (str | Path): Full path for the JSON output file.
    """
    with open(output_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    # [MATLAB] jsonencode(metrics) — single-line only; no indent support


def _save_run_report(metrics, optimizer_result, elapsed_sec, config, output_path):
    """Save a human-readable run report for cross-run comparison.

    Args:
        metrics (dict): Output of evaluate_metrics(...).
        optimizer_result (dict): Output of run_optimizer(...).
        elapsed_sec (float): Wall-clock optimization time in seconds.
        config (dict): Full parsed configuration dictionary.
        output_path (str | Path): Full path for the text output file.
    """
    converged     = optimizer_result["result"].success
    opt_cfg       = config.get("optimizer", {})
    directives    = config.get("directives", [])
    timestamp_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        "=" * 56,
        "Antenna Array Optimization Report",
        timestamp_str,
        "=" * 56,
        "",
        "RUN INFO",
        f"  Elapsed time:      {elapsed_sec:.1f} s",
        f"  Converged:         {converged}",
        f"  Iterations:        {metrics['iteration_count']}",
        f"  Final cost J:      {metrics['total_cost']:.6f}",
        "",
        "OPTIMIZER SETTINGS",
        "  Method:            L-BFGS-B",
        f"  Max iterations:    {opt_cfg.get('max_iterations')}",
        f"  Cost tolerance:    {opt_cfg.get('cost_tolerance')}",
        f"  N restarts:        {opt_cfg.get('n_restarts')}",
        f"  Amplitude bounds:  {opt_cfg.get('amplitude_bounds')}",
        f"  Phase only:        {opt_cfg.get('phase_only', False)}",
        f"  Amplitude only:    {opt_cfg.get('amplitude_only', False)}",
        "",
        "DIRECTIVES",
    ]
    for i, d in enumerate(directives):
        lines.append(
            f"  [{i}] {d['type']:<4}  theta={d['theta']:.1f}°  "
            f"phi={d.get('phi', 0.0):.1f}°  "
            f"width={d['width']:.1f}°  weight={d.get('weight', 1.0)}"
        )

    lines += [
        "",
        "RESULTS",
        f"  Global peak:       {metrics['global_peak_dbi']:.2f} dBi",
    ]
    if metrics["peak_to_null_ratio_db"] is not None:
        lines.append(f"  Peak-to-null:      {metrics['peak_to_null_ratio_db']:.2f} dB")
    lines.append("")

    for dm in metrics["directive_metrics"]:
        if dm["type"] == "peak":
            lines.append(
                f"  Peak @ theta={dm['theta_deg']:.1f}°, phi={dm['phi_deg']:.1f}°:  "
                f"{dm['gain_dbi']:.2f} dBi"
            )
        else:
            lines.append(
                f"  Null @ theta={dm['theta_deg']:.1f}°, phi={dm['phi_deg']:.1f}°:  "
                f"{dm['gain_dbi']:.2f} dBi  (depth = {dm['null_depth_db']:.2f} dB)"
            )
    lines += ["", "=" * 56]

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def _print_summary(metrics, optimizer_result):
    """Print a human-readable optimization summary to stdout.

    Args:
        metrics (dict): Output of evaluate_metrics(...).
        optimizer_result (dict): Output of run_optimizer(...).
    """
    converged = optimizer_result["result"].success
    print(f"  Converged:          {converged}")
    print(f"  Iterations:         {metrics['iteration_count']}")
    print(f"  Final cost J:       {metrics['total_cost']:.6f}")
    print(f"  Global peak:        {metrics['global_peak_dbi']:.2f} dBi")

    for dm in metrics["directive_metrics"]:
        if dm["type"] == "peak":
            print(f"  Peak @ theta={dm['theta_deg']}°, phi={dm['phi_deg']}°:  "
                  f"{dm['gain_dbi']:.2f} dBi")
        else:
            print(f"  Null @ theta={dm['theta_deg']}°, phi={dm['phi_deg']}°:  "
                  f"{dm['gain_dbi']:.2f} dBi  "
                  f"(depth = {dm['null_depth_db']:.2f} dB)")

    if metrics["peak_to_null_ratio_db"] is not None:
        print(f"  Peak-to-null:       {metrics['peak_to_null_ratio_db']:.2f} dB")


# ────────────────────────── MAIN ──────────────────────────────────

def main():
    """Run the full antenna array pattern optimization pipeline."""
    # ── Parse arguments ────────────────────────────────────────────
    arg_parser = argparse.ArgumentParser(
        description="Antenna Array Pattern Optimization"
    )
    arg_parser.add_argument(
        "--config", default="config.yaml",
        help="Path to the YAML configuration file (default: config.yaml).",
    )
    args = arg_parser.parse_args()

    # ── Load and validate config ───────────────────────────────────
    config = _load_config(args.config)
    _validate_config(config)

    # ── Create timestamped output directory ────────────────────────
    output_dir = _create_output_dir(config["output"]["results_dir"])
    print(f"Output directory: {output_dir}")

    # ── Stage 1: Load element patterns ────────────────────────────
    print("Loading element patterns...")
    patterns = load_element_patterns(config["element_patterns_dir"])

    # Pick the per-element complex pattern according to the chosen polarization.
    polarization = config.get("polarization", "copol")
    pattern_key  = POLARIZATION_TO_PATTERN_KEY[polarization]
    print(f"  Polarization: {polarization} (using '{pattern_key}')")

    # Stack individual element patterns into a single 3D array.
    element_patterns_stacked = np.stack([p[pattern_key] for p in patterns], axis=0)
    # [MATLAB] element_patterns_stacked = cat(3, patterns{:}.(pattern_key));
    theta_deg = patterns[0]["theta_deg"]
    phi_deg   = patterns[0]["phi_deg"]

    n_elements, n_theta, n_phi = element_patterns_stacked.shape
    print(f"  Loaded {n_elements} elements, grid {n_theta}×{n_phi}")

    # ── Stage 4: Run optimizer ─────────────────────────────────────
    print("Running optimizer...")
    directives    = config["directives"]
    t_start       = time.perf_counter()
    optimizer_result = run_optimizer(
        element_patterns_stacked, theta_deg, phi_deg,
        directives, config["optimizer"],
    )
    elapsed_sec     = time.perf_counter() - t_start
    weights_complex    = optimizer_result["weights_complex"]
    cost_history       = optimizer_result["cost_history"]
    all_cost_histories = optimizer_result["all_cost_histories"]
    best_run_index     = optimizer_result["best_run_index"]
    all_run_labels     = optimizer_result["all_run_labels"]

    # ── Stage 5: Evaluate metrics ──────────────────────────────────
    metrics = evaluate_metrics(
        element_patterns_stacked, theta_deg, phi_deg,
        weights_complex, directives, cost_history,
    )

    # ── Compute array factor dB grid for plotting ──────────────────
    # Array factor: coherent superposition of weighted element patterns.
    array_factor = compute_array_factor(weights_complex, element_patterns_stacked)
    power_linear = np.abs(array_factor) ** 2
    # Guard against log(0) before converting to dB.
    array_factor_db_grid = 10.0 * np.log10(np.maximum(power_linear, 1e-30))
    # [MATLAB] array_factor_db_grid = 10 * log10(max(abs(array_factor).^2, 1e-30));

    # ── Stage 6: Save plots ────────────────────────────────────────
    print("Saving plots...")
    save_all_plots(
        array_factor_db_grid, theta_deg, phi_deg,
        weights_complex, directives,
        all_cost_histories, best_run_index, all_run_labels,
        config["output"], output_dir,
    )

    # ── Stage 7: Save pattern evolution GIF (optional) ────────────
    if config["output"].get("save_pattern_gif", False):
        print("Saving pattern evolution GIF...")
        save_pattern_gif(
            optimizer_result["weights_history"],
            element_patterns_stacked,
            theta_deg, phi_deg, directives,
            cost_history, config["output"], output_dir,
        )

    # ── Save weights and metrics ───────────────────────────────────
    _save_weights_csv(weights_complex, output_dir / "weights.csv")
    _save_metrics_json(metrics,        output_dir / "metrics.json")

    # ── Save run report and config snapshot ───────────────────────
    _save_run_report(metrics, optimizer_result, elapsed_sec, config,
                     output_dir / "run_report.txt")
    shutil.copy(args.config, output_dir / "config.yaml")

    # ── Print summary ──────────────────────────────────────────────
    print("\nResults:")
    print(f"  Optimization time:  {elapsed_sec:.1f} s")
    _print_summary(metrics, optimizer_result)
    print(f"\nFiles saved to: {output_dir}")


if __name__ == "__main__":
    main()
