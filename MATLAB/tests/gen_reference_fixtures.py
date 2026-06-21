# ══════════════════════════════════════════════════════════════════
# GEN_REFERENCE_FIXTURES
# Generate Python reference outputs for the MATLAB port unit tests.
#
# Each `gen_<name>()` runs the existing src/ Python implementation on small,
# deterministic inputs and dumps the results to MATLAB/tests/fixtures/<name>.json.
# The matching MATLAB test (tests/test_<name>.m) loads the JSON and asserts the
# ported code reproduces these values.
#
# Usage:  python MATLAB/tests/gen_reference_fixtures.py [name1 name2 ...]
#         (no args = regenerate every fixture)
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════

import json
import sys
from pathlib import Path

import numpy as np

# Project root = three levels up from this file (MATLAB/tests/ -> repo root).
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
FIXTURES_DIR.mkdir(parents=True, exist_ok=True)


def _dump(name, payload):
    """Write payload to fixtures/<name>.json (lists/scalars only)."""
    out_path = FIXTURES_DIR / f"{name}.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    print(f"wrote {out_path}")


# ────────────────────────── STEP 1: weights_x ─────────────────────

def gen_weights_x():
    """Reference for x_to_weights / weights_to_x."""
    from src.cost.cost_function import x_to_weights, weights_to_x

    # A known complex weight vector with varied real/imag parts (incl. zeros, negatives).
    w = np.array([1.0 + 0.0j, 0.0 + 2.0j, -3.0 + 4.0j,
                  0.5 - 1.5j, -2.25 - 0.75j], dtype=complex)
    x_encoded = weights_to_x(w)

    # A separately chosen x vector to decode back into weights.
    x_in = np.array([0.0, 1.0, -2.0, 3.5, 4.0, -4.0], dtype=float)
    w_decoded = x_to_weights(x_in)

    _dump("weights_x", {
        "w_real":         w.real.tolist(),
        "w_imag":         w.imag.tolist(),
        "x":              x_encoded.tolist(),
        "x_in":           x_in.tolist(),
        "w_decode_real":  w_decoded.real.tolist(),
        "w_decode_imag":  w_decoded.imag.tolist(),
    })


# ────────────────────────── STEP 2: pure utilities ────────────────

def gen_nearest_index():
    """Reference for nearest_index (Python _nearest_index, 0-based argmin)."""
    from src.metrics.metrics import _nearest_index

    theta_deg = np.arange(0.0, 181.0, 5.0)            # 0..180 step 5
    targets   = [0.0, 2.4, 2.6, 7.5, 88.0, 180.0, 200.0, -10.0]
    idx0 = [int(_nearest_index(theta_deg, t)) for t in targets]   # 0-based

    _dump("nearest_index", {
        "grid_deg":     theta_deg.tolist(),
        "targets":      targets,
        "idx0_python":  idx0,    # MATLAB nearest_index should equal idx0 + 1
    })


def gen_angular_window_mask():
    """Reference for angular_window_mask on a small grid."""
    from src.cost.cost_function import angular_window_mask

    theta_deg = np.arange(0.0, 91.0, 10.0)     # 0..90 step 10  (10 pts)
    phi_deg   = np.arange(0.0, 360.0, 30.0)    # 0..330 step 30 (12 pts)

    cases = [
        {"tt": 30.0, "tp": 90.0, "tw": 20.0, "pw": 60.0},
        {"tt": 0.0,  "tp": 0.0,  "tw": 10.0, "pw": 10.0},
        {"tt": 50.0, "tp": 180.0, "tw": 25.0, "pw": 360.0},
    ]
    masks = []
    for c in cases:
        m = angular_window_mask(theta_deg, phi_deg, c["tt"], c["tp"], c["tw"], c["pw"])
        masks.append(m.astype(int).tolist())     # (N_theta x N_phi) of 0/1

    _dump("angular_window_mask", {
        "theta_deg": theta_deg.tolist(),
        "phi_deg":   phi_deg.tolist(),
        "cases":     cases,
        "masks":     masks,
    })


# ────────────────────────── STEP 3: array factor ──────────────────

def gen_compute_array_factor():
    """Reference for compute_array_factor on a synthetic complex stack."""
    from src.cost.cost_function import compute_array_factor

    rng = np.random.default_rng(42)
    n_elements, n_theta, n_phi = 3, 4, 5
    ep = (rng.standard_normal((n_elements, n_theta, n_phi))
          + 1j * rng.standard_normal((n_elements, n_theta, n_phi)))
    w = rng.standard_normal(n_elements) + 1j * rng.standard_normal(n_elements)

    af = compute_array_factor(w, ep)   # (n_theta, n_phi) complex

    _dump("compute_array_factor", {
        "shape":     [n_elements, n_theta, n_phi],
        "ep_real":   ep.real.tolist(),       # (N, N_theta, N_phi)
        "ep_imag":   ep.imag.tolist(),
        "w_real":    w.real.tolist(),
        "w_imag":    w.imag.tolist(),
        "af_real":   af.real.tolist(),       # (N_theta, N_phi)
        "af_imag":   af.imag.tolist(),
    })


# ────────────────────────── STEP 4: CST parser ────────────────────

# Directory of a real 4x4 CST export used for parser validation.
_CST_DIR = ROOT / "data" / "spacing0.9"


def gen_cst_parser():
    """Reference for parse_cst_file / load_element_patterns on real CST data."""
    from src.io.cst_parser import parse_cst_file, load_element_patterns, get_component

    # Parse one specific file (element [1]).
    one_file = _CST_DIR / "farfield (f=2400) [1].txt"
    p = parse_cst_file(one_file)
    copol = get_component(p, "Copol")
    cross = get_component(p, "Cross")

    n_theta = p["theta_deg"].size
    n_phi   = p["phi_deg"].size

    # Asymmetric sample points (0-based) that distinguish the theta vs phi axes.
    samples = [(0, 0), (1, 0), (0, 1), (5, 10), (90, 180),
               (n_theta - 1, n_phi - 1)]
    sample_vals = []
    for (t, q) in samples:
        sample_vals.append({
            "t": t, "q": q,
            "E_abs":         float(p["E_abs"][t, q]),
            "copol_abs":     float(copol["abs"][t, q]),
            "copol_phase":   float(copol["phase"][t, q]),
            "cross_abs":     float(cross["abs"][t, q]),
            "cross_phase":   float(cross["phase"][t, q]),
            "Ecx_real":      float(copol["complex"][t, q].real),
            "Ecx_imag":      float(copol["complex"][t, q].imag),
            "Xcx_real":      float(cross["complex"][t, q].real),
            "Xcx_imag":      float(cross["complex"][t, q].imag),
        })

    # Full slices: first phi-column (all theta) and first theta-row (all phi)
    # of the Copol complex field — pins down both axis orderings end-to-end.
    col0 = copol["complex"][:, 0]
    row0 = copol["complex"][0, :]

    # load_element_patterns: count + the parsed element ordering.
    patterns = load_element_patterns(_CST_DIR)

    _dump("cst_parser", {
        "n_theta":        int(n_theta),
        "n_phi":          int(n_phi),
        "theta_deg":      p["theta_deg"].tolist(),
        "phi_deg":        p["phi_deg"].tolist(),
        "samples":        sample_vals,
        "col0_real":      col0.real.tolist(),
        "col0_imag":      col0.imag.tolist(),
        "row0_real":      row0.real.tolist(),
        "row0_imag":      row0.imag.tolist(),
        "n_elements":     len(patterns),
        # Sanity scalar an element [1] copol_abs sum (order-independent checksum).
        "copol_abs_sum":  float(copol["abs"].sum()),
        "component_names": sorted(p["components"].keys()),
    })


# ────────────────────────── STEP 5-6: masks + cost ────────────────

# Shared small grid + directive set. The MATLAB tests replicate this directive
# list exactly (as a cell array of structs). Exercises normal windows,
# pole-crossing (theta<0 and theta>180), phi 0/360 wrap, and all aggregations.
def _small_grid():
    theta_deg = np.arange(0.0, 181.0, 10.0)    # 19 pts: 0..180
    phi_deg   = np.arange(0.0, 360.0, 30.0)    # 12 pts: 0..330
    return theta_deg, phi_deg


_DIRECTIVES = [
    {"type": "peak", "theta": 30.0,  "phi": 90.0,  "theta_width": 20.0, "phi_width": 60.0,  "weight": 1.0},
    {"type": "null", "theta": -30.0, "phi": 0.0,   "theta_width": 10.0, "phi_width": 20.0,  "weight": 3.0, "aggregation": "max"},
    {"type": "peak", "theta": 0.0,   "phi": 0.0,   "theta_width": 10.0, "phi_width": 10.0,  "weight": 2.0, "aggregation": "min"},
    {"type": "null", "theta": 90.0,  "phi": 350.0, "theta_width": 20.0, "phi_width": 40.0,  "weight": 5.0},
    {"type": "null", "theta": 185.0, "phi": 200.0, "theta_width": 20.0, "phi_width": 40.0,  "weight": 1.0},
]


def gen_phys_masks():
    """Reference physical masks for build_directive_physical_masks."""
    from src.cost.cost_function import build_directive_physical_masks

    theta_deg, phi_deg = _small_grid()
    masks = build_directive_physical_masks(theta_deg, phi_deg, _DIRECTIVES)
    masks_int = [m.astype(int).tolist() for m in masks]   # (n_dir, N_theta, N_phi)
    true_counts = [int(m.sum()) for m in masks]

    _dump("phys_masks", {
        "theta_deg":   theta_deg.tolist(),
        "phi_deg":     phi_deg.tolist(),
        "masks":       masks_int,
        "true_counts": true_counts,
    })


def gen_cost_function():
    """Reference J(x) values for build_cost_function across modes."""
    from src.cost.cost_function import build_cost_function

    theta_deg, phi_deg = _small_grid()
    n_theta, n_phi = theta_deg.size, phi_deg.size
    n_elements = 4

    rng = np.random.default_rng(7)
    ep  = (rng.standard_normal((n_elements, n_theta, n_phi))
           + 1j * rng.standard_normal((n_elements, n_theta, n_phi)))
    ep2 = (rng.standard_normal((n_elements, n_theta, n_phi))
           + 1j * rng.standard_normal((n_elements, n_theta, n_phi)))

    # Fixed variable vectors.
    x      = np.array([1.0, 0.2, -0.5, 0.7, 0.3, -0.9, 1.1, -0.4])  # 2N (standard/phase)
    x_amp  = np.array([0.8, 0.4, 1.2, 0.6])                          # N  (amplitude_only)

    j_standard = build_cost_function(ep, theta_deg, phi_deg, _DIRECTIVES, mode="standard")(x)
    j_phase    = build_cost_function(ep, theta_deg, phi_deg, _DIRECTIVES, mode="phase_only")(x)
    j_amp      = build_cost_function(ep, theta_deg, phi_deg, _DIRECTIVES, mode="amplitude_only")(x_amp)
    j_total    = build_cost_function(ep, theta_deg, phi_deg, _DIRECTIVES, mode="standard",
                                     element_patterns_secondary=ep2)(x)

    _dump("cost_function", {
        "theta_deg":  theta_deg.tolist(),
        "phi_deg":    phi_deg.tolist(),
        "ep_real":    ep.real.tolist(),
        "ep_imag":    ep.imag.tolist(),
        "ep2_real":   ep2.real.tolist(),
        "ep2_imag":   ep2.imag.tolist(),
        "x":          x.tolist(),
        "x_amp":      x_amp.tolist(),
        "j_standard": float(j_standard),
        "j_phase":    float(j_phase),
        "j_amp":      float(j_amp),
        "j_total":    float(j_total),
    })


# ────────────────────────── STEP 7: metrics ───────────────────────

_METRIC_DIRECTIVES = [
    {"type": "peak", "theta": 0.0,  "phi": 0.0, "theta_width": 20.0, "phi_width": 20.0, "weight": 1.0},
    {"type": "null", "theta": 40.0, "phi": 0.0, "theta_width": 10.0, "phi_width": 20.0, "weight": 2.0},
]


def gen_hpbw():
    """Reference 3 dB beamwidths (grid-index units) for _compute_hpbw."""
    from src.metrics.metrics import _compute_hpbw

    # (a) symmetric peak, no pole wrap.
    a_cut   = np.array([-20., -10., -4., -2., 0., -2., -4., -10., -20.])
    a_peak0 = 4
    a_width = _compute_hpbw(a_cut, a_peak0)

    # (b) peak at the left edge -> left side wraps into the opposite cut.
    b_cut   = np.array([0., -2., -6., -12., -20., -25., -28., -30.])
    b_opp   = np.array([-1., -4., -10., -18., -25., -30., -33., -35.])
    b_peak0 = 0
    b_width = _compute_hpbw(b_cut, b_peak0, opposite_cut=b_opp)

    # (c) peak at the right edge -> right side wraps into the opposite cut.
    c_cut   = np.array([-30., -28., -25., -20., -12., -6., -2., 0.])
    c_opp   = np.array([-35., -33., -30., -25., -18., -10., -4., -1.])
    c_peak0 = len(c_cut) - 1
    c_width = _compute_hpbw(c_cut, c_peak0, opposite_cut=c_opp)

    _dump("hpbw", {
        "a_cut": a_cut.tolist(), "a_peak0": a_peak0, "a_width": float(a_width),
        "b_cut": b_cut.tolist(), "b_opp": b_opp.tolist(), "b_peak0": b_peak0, "b_width": float(b_width),
        "c_cut": c_cut.tolist(), "c_opp": c_opp.tolist(), "c_peak0": c_peak0, "c_width": float(c_width),
    })


def gen_directivity_grid():
    """Reference directivity dBi grid for _compute_directivity_dbi_grid."""
    from src.metrics.metrics import _compute_directivity_dbi_grid

    theta_deg = np.arange(0.0, 121.0, 30.0)    # 5 pts
    phi_deg   = np.arange(0.0, 360.0, 60.0)    # 6 pts
    rng = np.random.default_rng(11)
    af = (rng.standard_normal((5, 6)) + 1j * rng.standard_normal((5, 6)))

    dbi_auto = _compute_directivity_dbi_grid(af, theta_deg, phi_deg)
    dbi_norm = _compute_directivity_dbi_grid(af, theta_deg, phi_deg, normalizer_power=2.5)

    _dump("directivity_grid", {
        "theta_deg": theta_deg.tolist(),
        "phi_deg":   phi_deg.tolist(),
        "af_real":   af.real.tolist(),
        "af_imag":   af.imag.tolist(),
        "dbi_auto":  dbi_auto.tolist(),
        "dbi_norm":  dbi_norm.tolist(),
    })


def gen_evaluate_metrics():
    """Reference scalar metrics on real CST data (uniform weights, 16-elem URA)."""
    from src.io.cst_parser import load_element_patterns, get_component
    from src.metrics.metrics import evaluate_metrics

    patterns = load_element_patterns(_CST_DIR)
    theta_deg = patterns[0]["theta_deg"]
    phi_deg   = patterns[0]["phi_deg"]
    ep = np.stack([get_component(p, "Copol")["complex"] for p in patterns], axis=0)
    n = ep.shape[0]
    w = np.ones(n, dtype=complex)   # uniform -> broadside beam (tests pole wrap HPBW)

    m = evaluate_metrics(ep, theta_deg, phi_deg, w, _METRIC_DIRECTIVES, [])

    dms = [{"type": d["type"], "gain_dbi": d["gain_dbi"],
            "null_depth_db": d["null_depth_db"], "cost_term": d["cost_term"]}
           for d in m["directive_metrics"]]

    _dump("evaluate_metrics", {
        "n_elements":            n,
        "global_peak_dbi":       m["global_peak_dbi"],
        "global_peak_theta_deg": m["global_peak_theta_deg"],
        "global_peak_phi_deg":   m["global_peak_phi_deg"],
        "hpbw_theta_deg":        m["hpbw_theta_deg"],
        "hpbw_phi_deg":          m["hpbw_phi_deg"],
        "total_cost":            m["total_cost"],
        "peak_to_null_ratio_db": m["peak_to_null_ratio_db"],
        "directive_metrics":     dms,
    })


# ────────────────────────── STEP 8: windows ───────────────────────

def gen_windows():
    """Reference values for all named 1-D windows at several lengths."""
    from scripts.compare_classical import _window_1d, CLASSICAL_TECHNIQUE_NAMES

    cases = []
    for n in (3, 4, 5, 8):
        windows = {name: _window_1d(name, n).tolist() for name in CLASSICAL_TECHNIQUE_NAMES}
        cases.append({"n": n, "windows": windows})
    _dump("windows", {"names": list(CLASSICAL_TECHNIQUE_NAMES), "cases": cases})


# ────────────────────────── STEP 9: optimizer ─────────────────────

def gen_optimizer():
    """Reference converged cost for run_optimizer on a small synthetic URA."""
    from scripts.compare_classical import build_ura_element_patterns
    from src.optimize.optimizer import run_optimizer
    from src.cost.cost_function import build_cost_function, weights_to_x

    n_side = 3
    theta_deg = np.arange(0.0, 181.0, 10.0)
    phi_deg   = np.arange(0.0, 360.0, 30.0)
    ep = build_ura_element_patterns(n_side, 0.5, theta_deg, phi_deg)   # (9, 19, 12)
    n_elements = ep.shape[0]

    directive = {"type": "peak", "theta": 20.0, "phi": 60.0,
                 "theta_width": 30.0, "phi_width": 60.0, "weight": 1.0}
    directives = [directive]

    config = {
        "max_iterations": 200,
        "cost_tolerance": 1e-9,
        "n_restarts": 2,
        "amplitude_bounds": [0.0, 1.0],
        "use_single_element_init": True,
    }

    opt = run_optimizer(ep, theta_deg, phi_deg, directives, config)
    j_star  = float(opt["result"].fun)
    n_runs  = len(opt["all_cost_histories"])

    # Cost at uniform weights (the optimizer must do at least this well).
    cost_fn = build_cost_function(ep, theta_deg, phi_deg, directives, mode="standard")
    j_uniform = float(cost_fn(weights_to_x(np.ones(n_elements, dtype=complex))))

    _dump("optimizer", {
        "theta_deg":  theta_deg.tolist(),
        "phi_deg":    phi_deg.tolist(),
        "ep_real":    ep.real.tolist(),
        "ep_imag":    ep.imag.tolist(),
        "n_elements": n_elements,
        "directive":  directive,
        "config":     config,
        "j_star":     j_star,
        "j_uniform":  j_uniform,
        "n_runs":     n_runs,
    })


# ────────────────────────── STEP 10: YAML config ──────────────────

def gen_config_yaml():
    """Reference values from PyYAML parsing of the two real config files."""
    import yaml

    with open(ROOT / "config.yaml", "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    with open(ROOT / "test_config.yaml", "r", encoding="utf-8") as fh:
        tcfg = yaml.safe_load(fh)

    config_checks = {
        "element_patterns_dir": cfg["element_patterns_dir"],
        "polarization":         cfg["polarization"],
        "n_directives":         len(cfg["directives"]),
        "d0_type":              cfg["directives"][0]["type"],
        "d0_theta":             cfg["directives"][0]["theta"],
        "d0_phi":               cfg["directives"][0]["phi"],
        "d0_theta_width":       cfg["directives"][0]["theta_width"],
        "d0_phi_width":         cfg["directives"][0]["phi_width"],
        "d0_weight":            cfg["directives"][0]["weight"],
        "d0_level_is_null":     cfg["directives"][0]["level"] is None,
        "d1_type":              cfg["directives"][1]["type"],
        "d1_theta":             cfg["directives"][1]["theta"],
        "d1_weight":            cfg["directives"][1]["weight"],
        "opt_max_iter":         cfg["optimizer"]["max_iterations"],
        "opt_ctol":             cfg["optimizer"]["cost_tolerance"],
        "opt_n_restarts":       cfg["optimizer"]["n_restarts"],
        "opt_bounds":           cfg["optimizer"]["amplitude_bounds"],
        "opt_phase_only":       cfg["optimizer"]["phase_only"],
        "opt_amp_only":         cfg["optimizer"]["amplitude_only"],
        "opt_use_uniform":      cfg["optimizer"]["use_uniform_init"],
        "opt_use_single":       cfg["optimizer"]["use_single_element_init"],
        "out_results_dir":      cfg["output"]["results_dir"],
        "out_cut_type":         cfg["output"]["plot_cut_type"],
        "out_dyn_range":        cfg["output"]["plot_dynamic_range_db"],
        "out_save_gif":         cfg["output"]["save_pattern_gif"],
        "out_gif_frames":       cfg["output"]["gif_max_frames"],
        "out_equal_area":       cfg["output"]["plot_equal_area"],
    }

    sc0 = tcfg["scenarios"][0]
    test_checks = {
        "n_side":             tcfg["n_side"],
        "d_over_lambda":      tcfg["d_over_lambda"],
        "theta_step_deg":     tcfg["theta_step_deg"],
        "phi_step_deg":       tcfg["phi_step_deg"],
        "element_source":     tcfg["element_source"],
        "element_patterns_dir": tcfg["element_patterns_dir"],
        "polarization":       tcfg["polarization"],
        "n_restarts":         tcfg["n_restarts"],
        "max_iterations":     tcfg["max_iterations"],
        "cost_tolerance":     tcfg["cost_tolerance"],
        "plot_dynamic_range_db": tcfg["plot_dynamic_range_db"],
        "n_scenarios":        len(tcfg["scenarios"]),
        "s0_label":           sc0["label"],
        "s0_steer_theta":     sc0["steer_theta_deg"],
        "s0_null_theta":      sc0["null_theta_deg"],
        "s0_n_directives":    len(sc0["directives"]),
        "s0_d0_type":         sc0["directives"][0]["type"],
        "s0_d0_theta":        sc0["directives"][0]["theta"],
        "s0_d0_phi":          sc0["directives"][0]["phi"],
        "s0_d0_phi_width":    sc0["directives"][0]["phi_width"],
        "s0_d1_type":         sc0["directives"][1]["type"],
        "s0_d1_theta":        sc0["directives"][1]["theta"],
        "s0_d1_weight":       sc0["directives"][1]["weight"],
    }

    _dump("config_yaml", {"config": config_checks, "test_config": test_checks})


# ────────────────────────── STEP 11: compare_classical helpers ────

def gen_compare_helpers():
    """Reference for the deterministic compare_classical math helpers."""
    from scripts.compare_classical import (
        build_ura_element_patterns, _steering_phase_vector,
        _data_driven_steering_vector, classical_weights,
        principal_plane_cut, principal_plane_theta_axis,
        _power_normalize_weights, evaluate_directive_metrics,
    )
    from src.cost.cost_function import compute_array_factor

    n_side, d = 3, 0.5
    theta_deg = np.arange(0.0, 181.0, 10.0)
    phi_deg   = np.arange(0.0, 360.0, 30.0)
    ep = build_ura_element_patterns(n_side, d, theta_deg, phi_deg)

    st, sp = 20.0, 60.0
    steering_geo = _steering_phase_vector(n_side, d, st, sp)
    steering_dd  = _data_driven_steering_vector(ep, theta_deg, phi_deg, st, sp)

    names = ["uniform", "hamming", "taylor_25dB", "chebyshev_25dB"]
    cw = {}
    for nm in names:
        w = classical_weights(nm, n_side, d, st, sp, ep, theta_deg, phi_deg)
        cw[nm] = {"real": w.real.tolist(), "imag": w.imag.tolist()}

    # Cut + directive metrics on a concrete pattern (hamming, normalised).
    w_ham   = classical_weights("hamming", n_side, d, st, sp, ep, theta_deg, phi_deg)
    w_norm  = _power_normalize_weights(w_ham)
    af      = compute_array_factor(w_norm, ep)
    af_mag  = np.abs(af)
    cut     = principal_plane_cut(af, phi_deg)
    directive = {"type": "peak", "theta": 20.0, "phi": 60.0,
                 "theta_width": 30.0, "phi_width": 60.0, "weight": 1.0}
    edm = evaluate_directive_metrics(af_mag, theta_deg, phi_deg, directive)

    # power_normalize standalone.
    w_test = np.array([1 + 1j, 2 - 0.5j, -3 + 0j, 0 + 4j], dtype=complex)
    w_test_norm = _power_normalize_weights(w_test)

    _dump("compare_helpers", {
        "theta_deg": theta_deg.tolist(),
        "phi_deg":   phi_deg.tolist(),
        "ep_real":   ep.real.tolist(),
        "ep_imag":   ep.imag.tolist(),
        "steer_theta": st, "steer_phi": sp,
        "steering_geo_real": steering_geo.real.tolist(),
        "steering_geo_imag": steering_geo.imag.tolist(),
        "steering_dd_real":  steering_dd.real.tolist(),
        "steering_dd_imag":  steering_dd.imag.tolist(),
        "cw_names": names,
        "cw": cw,
        "theta_axis": principal_plane_theta_axis(theta_deg).tolist(),
        "af_mag":    af_mag.tolist(),
        "cut":       cut.tolist(),
        "directive": directive,
        "edm": {"peak_dbi": edm["peak_dbi"], "min_db": edm["min_db"],
                "mean_db": edm["mean_db"], "max_db": edm["max_db"]},
        "wtest_real": w_test.real.tolist(), "wtest_imag": w_test.imag.tolist(),
        "wtest_norm_real": w_test_norm.real.tolist(),
        "wtest_norm_imag": w_test_norm.imag.tolist(),
    })


def gen_run_scenario():
    """Reference classical-technique results from run_scenario (deterministic)."""
    from scripts.compare_classical import (
        build_ura_element_patterns, run_scenario, CLASSICAL_TECHNIQUE_NAMES,
    )

    n_side, d = 3, 0.5
    theta_deg = np.arange(0.0, 181.0, 10.0)
    phi_deg   = np.arange(0.0, 360.0, 30.0)
    ep = build_ura_element_patterns(n_side, d, theta_deg, phi_deg)

    scenario = {
        "steer_theta_deg": 0.0, "steer_phi_deg": 0.0,
        "directives": [
            {"type": "peak", "theta": 0.0,  "phi": 0.0, "theta_width": 20.0, "phi_width": 360.0, "weight": 1.0},
            {"type": "null", "theta": 40.0, "phi": 0.0, "theta_width": 10.0, "phi_width": 360.0, "weight": 2.0},
        ],
    }
    cfg = {"n_restarts": 2, "max_iterations": 100, "cost_tolerance": 1e-9}
    results = run_scenario(scenario, ep, theta_deg, phi_deg, n_side, d, cfg)

    classical = {}
    for nm in CLASSICAL_TECHNIQUE_NAMES:
        r = results[nm]
        classical[nm] = {
            "cut": r["cut_vm"].tolist(),
            "dm": [{"peak_dbi": dm["peak_dbi"], "min_db": dm["min_db"],
                    "mean_db": dm["mean_db"], "max_db": dm["max_db"]}
                   for dm in r["directive_metrics"]],
        }

    _dump("run_scenario", {
        "names": list(CLASSICAL_TECHNIQUE_NAMES),
        "n_cut": len(results["uniform"]["cut_vm"]),
        "n_directives": len(scenario["directives"]),
        "classical": classical,
        "has_optimized": "optimized" in results,
    })


# ────────────────────────── DISPATCH ──────────────────────────────

GENERATORS = {
    "weights_x": gen_weights_x,
    "nearest_index": gen_nearest_index,
    "angular_window_mask": gen_angular_window_mask,
    "compute_array_factor": gen_compute_array_factor,
    "cst_parser": gen_cst_parser,
    "phys_masks": gen_phys_masks,
    "cost_function": gen_cost_function,
    "hpbw": gen_hpbw,
    "directivity_grid": gen_directivity_grid,
    "evaluate_metrics": gen_evaluate_metrics,
    "windows": gen_windows,
    "optimizer": gen_optimizer,
    "config_yaml": gen_config_yaml,
    "compare_helpers": gen_compare_helpers,
    "run_scenario": gen_run_scenario,
}


def main():
    requested = sys.argv[1:] or list(GENERATORS.keys())
    for name in requested:
        if name not in GENERATORS:
            raise SystemExit(f"Unknown fixture '{name}'. Known: {list(GENERATORS)}")
        GENERATORS[name]()


if __name__ == "__main__":
    main()
