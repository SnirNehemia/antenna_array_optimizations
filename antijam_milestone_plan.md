# Milestone plan — Adaptive anti-jam beamforming (Option B: classical core + non-stationary bandit)

**Status:** ACTIVE — implementation started 2026-07-18 (P0)
**Owner:** Snir
**Cadence assumption:** 5–10 hrs/week, Claude Code assisted
**Target duration:** 8–9 weeks calendar
**Implementation language:** **MATLAB only**, directly on top of the existing `MATLAB/` port
(R2020a + Optimization Toolbox — see Section 8). The Python side stays frozen at Milestone 1.
**Living document:** update the `Status` field of each phase (`todo / in-progress / done / blocked`) as work proceeds. Append session summaries to `docs/notes.md` per project convention; reference the phase ID (P0–P6) in each entry.

---

## 1. Objective

Extend the existing pattern-synthesis tool (Milestone 1: user-specified peaks/nulls → optimized complex weights) into a **closed-loop adaptive system** that maintains a radiation peak toward a predefined desired direction while automatically nulling a **single, unknown, slowly moving jammer**, using only runtime measurements. Implemented **directly in MATLAB** in a new `MATLAB/antijam_utils/` function library, reusing the validated Milestone-1 port in `MATLAB/matlab_utils/` (parser, cost function, fmincon optimizer, YAML reader, plotting helpers) without modifying it.

### Locked scope decisions
| Item | Decision |
|---|---|
| Language | MATLAB only (no Python counterpart; no migration phase) |
| Jammers | 1 simultaneous, unknown angle, slow drift (seconds timescale), power steps / on-off allowed |
| Geometry | ~~2D azimuth cut~~ **AMENDED [P7], 2026-07-20: native 2-D (theta, phi) tracking** — see P7 below |
| Element data | CST embedded element patterns via `load_element_patterns` + `stack_component` → `(N_el, N_theta, N_phi)` complex stack; the anti-jam engine works on a principal-plane cut → `E` matrix `(N_el × N_angles)` complex. Common phase center; no geometric phase re-added. |
| Simulation fidelity | Pattern-level: SINR computed from pattern gains + jammer/noise powers. No waveform/IQ modeling. |
| Observation modes | **Mode C (covariance):** per-element snapshots available → classical LCMV/MVDR. **Mode S (scalar):** only output SINR feedback → bandit + SPSA. |
| RL approach | Non-stationary multi-armed bandit over a codebook of pre-synthesized weight vectors. No neural networks. |
| Deliverables | `MATLAB/antijam_utils/` library + closed-loop simulation + KPI evaluation + automated report script |

### Algorithm inputs (per customer spec)
1. Desired peak direction θ_s (deg) and operating SINR threshold `sinr_min_db` ("system still operates").
2. Measured element patterns: CST export, same format/parser as Milestone 1 (`parse_cst_file` / `load_element_patterns` reused as-is).

---

## 2. Core physics (reference formulas)

With `e(θ)` = column of `E` at angle θ (complex embedded element responses; steering vector — mutual coupling included by construction):

- **Output SINR:**
  `SINR(w) = σ_s² |wᴴ e(θ_s)|² / ( σ_j² |wᴴ e(θ_j)|² + σ_n² ‖w‖² )`
  (per-element noise assumed i.i.d. → noise power scales with ‖w‖².)
- **Interference-plus-noise covariance (analytic):**
  `R = σ_j² e(θ_j) e(θ_j)ᴴ + σ_n² I`
- **MVDR/LCMV oracle solution:**
  `w_opt = R⁻¹ e(θ_s) / ( e(θ_s)ᴴ R⁻¹ e(θ_s) )`
  → defines the **oracle upper bound** for all algorithms.
- **Snapshot model (Mode C measurement):**
  `x(t) = e(θ_s) s(t) + e(θ_j) j(t) + n(t)`, with `s, j, n` circular complex Gaussian
  (sampled as `(randn + 1i*randn)/sqrt(2)` scaled by the component's σ).
  Sample covariance with exponential forgetting: `R̂(t) = λ R̂(t−1) + (1−λ) x(t) x(t)ᴴ`.

All optimization in linear field magnitude; display in dB. Phase radians internally, degrees for display (project convention).

---

## 3. Module breakdown

New code in **bold**; Milestone-1 `matlab_utils/` functions reused **unmodified** (any needed change to them must be noted here first, per `CLAUDE.md`). Layout follows the repo convention: flat folder, file = function, module boundary expressed by filename prefix (`sim_` / `adapt_` / `agent_` / `kpi_` / `plot_`).

```
MATLAB/
├── matlab_utils/               # REUSED, unmodified:
│                               #   parse_cst_file, load_element_patterns, stack_component,
│                               #   principal_plane_cut, run_optimizer, build_cost_function,
│                               #   read_config_yaml, nearest_index, power_normalize_weights, ...
├── **antijam_utils/**          # NEW — flat, file = function
│   ├── sim_scenario.m          # jammer trajectory generator: θ_j(t) drift, power steps, on/off
│   ├── sim_engine_init.m       # precompute E cut, steering columns, powers → sim state struct
│   ├── sim_engine_step.m       # closed-loop stepper (the ONLY channel algorithms see)
│   ├── sim_analytic_covariance.m  # R(θ_j) for the oracle
│   ├── adapt_lcmv.m            # MVDR/LCMV weights from R̂ + steering vector; diagonal loading
│   ├── adapt_tracking_init.m   # ┐ Mode C tracker: recursive covariance estimation
│   ├── adapt_tracking_update.m # ┘ (exponential forgetting) + LCMV weight recompute
│   ├── adapt_spsa_init.m       # ┐ Mode S baseline: gradient-free SPSA
│   ├── adapt_spsa_update.m     # ┘ (2 SINR probes per update, warm-startable)
│   ├── agent_codebook_build.m  # arm generation via run_optimizer (peak θ_s + null sector θ_k)
│   ├── agent_bandit_init.m     # ┐ discounted Thompson sampling (primary),
│   ├── agent_bandit_update.m   # ┘ sliding-window UCB (config alternative)
│   ├── kpi_evaluate.m          # the 5 milestone KPIs (Section 5) from timeline logs
│   ├── plot_antijam_report.m   # SINR timeline, arm-selection heatmap, pattern snapshots
│   ├── select_polarization_stacks.m  # ┐ shared plumbing (run_antijam + demo):
│   ├── extract_cut.m                 # │ polarization stacks, 1-D cut,
│   ├── closed_loop_run.m             # │ one algorithm × scenario closed loop,
│   ├── write_kpi_table.m             # ┘ KPI aggregation to CSV/txt
│   ├── save_run_gif.m          # animated 2-D pattern + jammer dot + live SINR/gain traces
│   ├── run_jammer_demo.m       # jammer_config.yaml-driven demo: scenarios × algorithms → GIFs
│   └── run_antijam.m           # closed-loop campaign driver: scenario × algorithm × KPI report
├── scripts/
│   ├── run_antijam_script.m    # thin entry-point wrappers (addpath matlab_utils
│   └── run_jammer_demo_script.m #  + antijam_utils)
└── tests/
    └── test_antijam_*.m        # analytic / toy-case unit tests (Section 8, Validation)
```

`config.yaml` is EXTENDED with `antijam / sim / adapt / agent` sections (schema in Section 6), parsed by the existing `read_config_yaml`.

**Interface contracts (freeze in P0):**
- `[obs, sim_state] = sim_engine_step(sim_state, w)` — `obs` is a struct with fields `sinr_db` (always) and `snapshots` (`N_el × K` complex; Mode C only, `[]` in Mode S). This is the *only* channel algorithms see; the Mode C / Mode S distinction is enforced by which `obs` fields the algorithm is allowed to consume. θ_j never leaves the sim state.
- All algorithms expose a common signature `[w, alg_state] = <alg>_update(alg_state, obs)` so the runner and KPI code are algorithm-agnostic. State is a plain struct created by a matching `<alg>_init(config, ...)` function.
- `W = agent_codebook_build(theta_s_deg, null_grid_deg, config)  % (N_el × n_arms)` — pure reuse of `run_optimizer` with peak + null-sector directives; result cached to a `.mat` file and reproducible from `config.yaml`.

---

## 4. Phased plan with success criteria

### P0 — Interfaces & spec freeze (week 1) — Status: **done** (2026-07-18: signatures approved in-session; config schema parse-verified)
Define the contracts above, extend the `config.yaml` schema, add pointer to this plan in `CLAUDE.md`, stub all new functions with header comments only (repo docstring style). Verify `read_config_yaml` parses the extended schema (nested block mappings — see Section 6 note).
**Success criteria:** interface signatures reviewed and approved in-session before any implementation; open decisions in Section 7 resolved or explicitly deferred with defaults.

### P1 — Simulation harness (weeks 1–3) — Status: **done** (2026-07-18: all four gates pass — null depth −75.9 dB, cov error 1.3% @ N=400, drift 0.07% error, toy SINR 1e-10; `tests/test_antijam_sim.m`. `adapt_lcmv` implemented early, needed by the oracle gate.)
Implement `sim_scenario`, `sim_engine_init`, `sim_engine_step`, `sim_analytic_covariance` + unit tests.
**Success criteria (quantitative):**
- Oracle LCMV null depth at θ_j ≥ 40 dB below peak for σ_j/σ_n = 20 dB (analytic-R case), matching theory.
- Sample covariance from snapshots converges to analytic R: relative Frobenius error < 5% at N = 50·N_el snapshots.
- Scenario generator reproduces configured drift rate within ±2% over a 60 s trajectory; deterministic under fixed seed (`rng(seed, 'twister')`).
- SINR from the engine matches a hand-computed value on a 2-element toy case to 1e-10 relative.

### P2 — Mode C: covariance adaptation (weeks 3–4.5) — Status: **done** (2026-07-18: λ/loading swept — chose λ=0.98, loading **+10 dB** (MPDR snapshots self-null under light loading; −10 dB default gave 2–3.7 dB gap). Gates: gap 0.39–0.80 dB across S1–S5 (< 1), recovery 1 update (≤ 25), availability 97.4% on S2 (≥ 95%). `tests/test_antijam_tracking.m`.)
`adapt_lcmv` + `adapt_tracking_update`; diagonal loading tuned; forgetting factor λ swept.
**Success criteria:**
- Steady-state oracle gap < 1 dB on the standard scenario suite (Section 5).
- Recovery time after a 10° jammer jump ≤ 25 covariance updates.
- SINR availability ≥ 95% (threshold per config) against the drifting-jammer scenario.

### P3 — Mode S baseline: SPSA (week 5) — Status: **done** (2026-07-18: gate passed — median 10 probes to oracle−3 dB, 50/50 seeds (≤ 150 median). Tuning finding: stability constant `A=15` + step-norm clamp `step_max=0.3` are load-bearing (dB-scale gradients explode near nulls); both added as required `adapt.spsa` keys. Drift characterization saved to `results/antijam/p3_spsa_drift_characterization.png` — availability 34–64% under drift vs 97% for Mode C, motivating P5. `tests/test_antijam_spsa.m`.)
`adapt_spsa_update`, warm-startable from any weight vector.
**Success criteria:**
- Static jammer: reaches within 3 dB of oracle SINR in ≤ 150 SINR probes (median over 50 Monte Carlo seeds).
- Documented characterization (plot + note in `docs/notes.md`) of tracking behavior vs drift rate — expected to degrade; this motivates the bandit.

### P4 — Codebook generation (weeks 5–6) — Status: **done** (2026-07-19: gates pass on a 16-el ULA toy at real-array proportions; `tests/test_antijam_codebook.m`; coverage plot `results/antijam/p4_codebook_coverage.png`. **Design findings:** (1) peak directive needs `aggregation: 'min'` with width ≈ HPBW, else the solid-angle-weighted mean parks the beam off θ_s; (2) the composite cost plateaus at ~−25 dB null depth (true Pareto point) — arms are polished by **null-space projection** onto the window steering columns → numerically exact nulls at sampled angles; per-arm gate refined to depth at the null CENTER ≤ −30 dB + the −25 dB coverage gate (window-MEAN depth is DOF-limited for filled nulls and is reported, not gated); (3) `guard_deg` must exceed ≈ (peak_width+null_width)/2 AND keep windows ~4+ beamwidths from boresight or the projection eats the main lobe — builder warns on the former.)
`agent_codebook_build` calling `run_optimizer` once per candidate null sector (grid resolution per config, default = CST angular step or coarser).
**Success criteria:**
- Every arm: null-sector depth ≤ −30 dB relative to peak AND peak gain within 1 dB of the unconstrained optimum at θ_s.
- Coverage plot generated: worst-case achievable null depth vs jammer angle over the full azimuth span (no coverage hole worse than −25 dB between adjacent arms).
- Codebook build is offline, cached to a `.mat` file, and reproducible from `config.yaml`.

### P5 — Bandit agent (weeks 6–8) — Status: **done** (2026-07-19: gates pass — static identification 93/100 in ≤ 30 probes, S2 availability 95.1% median (≥ 90%), recovery median 0 probes vs SPSA's 25 (≥ 2×); `tests/test_antijam_bandit.m`. **Gate-1 refinement:** arms whose natural sidelobe nulls coincide with θ_j earn rewards identical to the designated arm, so success = window coverage OR within 1 dB of the best arm's SINR. **Tuning:** posterior scale σ̃ = 1 dB is load-bearing (5 dB → 40% identification); re-exploration under drift comes from the discount, not σ̃. swucb alternative implemented + smoke-tested. Stretch goal (SPSA fine-tune within arm) not needed — deferred per descope order. Regret curve deferred to the P6 report.)
`agent_bandit_update`: **discounted Thompson sampling** with Gaussian reward model on measured SINR (primary); sliding-window UCB as a config-selectable alternative. Optional: SPSA fine-tune within the winning arm's sector (stretch goal, cut first if schedule slips).
**Success criteria:**
- Static jammer: correct arm (or a neighbor covering the true angle) selected within ≤ 30 probes with ≥ 90% success over 100 Monte Carlo runs.
- Drifting jammer: SINR availability ≥ 90% on the standard suite; arm-selection heatmap visibly tracks θ_j(t).
- Beats SPSA baseline on recovery time after jammer jumps by ≥ 2× (median).
- Regret curve plotted vs oracle for the report.

### P6 — Evaluation campaign & demo (weeks 8–9) — Status: **done** (2026-07-19, gates partially met — see below)
Monte Carlo runner over the scenario suite × {LCMV, SPSA, bandit, oracle}; full KPI table; `run_antijam` + `run_antijam_script` producing the headline plots.
**Success criteria:**
- Automated report (single command) writing to a timestamped `results/antijam/<ts>/` folder: KPI table + plots — SINR timeline with threshold band, recovery-time histograms, arm heatmap, oracle-gap bars, pattern snapshots at selected instants. — ✔ **met** (24 artifacts incl. regret curves and the S6 null-lifecycle figure; `kpi_table.csv/txt`).
- Headline metric met: SINR availability ≥ 90% (Mode S) / ≥ 95% (Mode C) at the configured `sinr_min_db` across the suite. — **Mode C ✔ met** (LCMV 97.4–99.9% across S1–S6; oracle gap 0.6–1.5 dB). **Mode S partially met**: bandit ≥ 90% on S1/S4/S5/S6 (92.2/92.1/90.1/97.2%) but **85–86% on the sustained-drift cases S2/S3** — diagnosed as irreducible exploration overhead of the vanilla discounted-TS at the 5 dB operating headroom (85% of below-threshold steps are exploration probes of non-covering arms while a ≥ 11 dB arm always exists). Forward paths (not in scope): neighbor-restricted exploration, higher link headroom, or LCMV/bandit hybrid. SPSA baseline is slow on the real 20-element array (~1000 probes to converge even warm-started from the quiescent beam) — retained as the documented motivation for the bandit, per P3.

**P6 real-data findings (ManyDipoles, Theta pol, φ-cut @ θ=90°):** quiescent gain 12.3 dB, HPBW ≈ 15° → `sigma_s_db: 3` (15.3 dB jammer-free SINR), `guard_deg: 45`. Codebook tuning: 20°-wide windows on a 10° grid (28 arms) so the arm count stays inside the discount horizon (1/(1−0.99) = 100) and drift handoffs halve; the **null-space projection is restricted to the cut's own steering columns** — projecting the full 2-D θ×φ window ate rank 15 of 20 DOF and destroyed the beam. Bandit `discount: 0.99`, `sigma_tilde_db: 1` (new required key).

**Schedule buffer:** ~1 week implicit slack. **Descope order if needed:** P5 stretch goal (SPSA fine-tune) → sliding-window UCB alternative → reduce Monte Carlo scenario count in P6.

### P7 — Native 2D (theta, phi) engine (amendment) — Status: **done** (2026-07-20, demo-scoped; see below)
Triggered by a real request the 1-D cut couldn't express: target and jammer at
independently different theta AND phi. Full plan:
`C:\Users\snirn\.claude\plans\giggly-wishing-whisper.md`; session-log detail in
`docs/notes.md` `[P7]` entry (2026-07-20).

Replaced the 1-D cut end-to-end: `sim_scenario`/`sim_engine_init/step`/
`sim_analytic_covariance`/`closed_loop_run` operate on the full `(theta, phi)`
grid natively; `agent_codebook_build` synthesizes nulls on a genuine 2-D
candidate grid with a rank-capped projection (new `agent.null_rank_cap` key,
replacing the old cut-restriction trick); `kpi_evaluate`'s null-pointing KPI
and the report/GIF jammer-dot plotting are natively 2-D. `extract_cut.m`
deleted. Config schema: `theta_s_deg` + `phi_s_deg` replace `cut_type` /
`cut_theta_deg` / `cut_phi_deg`; scenarios take `theta_j_deg` + `phi_j_deg` (or
per-axis drift rates) instead of a single `angle_deg`.

**What's validated:** all `test_antijam_*.m` unit/gate tests pass under the
new 2-D contracts (some toy-array gate thresholds relaxed for 2-D codebook
density — see `docs/notes.md`). `run_jammer_demo` end-to-end on real CST data
(`data/patchs_with_monopoles`, target θ=30°/φ=90°, jammer θ=45°/φ=180°) runs
to completion and renders the jammer at its true 2-D position.

**What's NOT validated (follow-up, not blocking):** the full P1–P6 Monte
Carlo campaign (`run_antijam` on `config.yaml`'s real-array scenario suite)
has not been rerun under 2-D geometry — the quantitative KPI numbers quoted
in P1–P6 above are from the 1-D-cut engine and are stale for a 2-D reading.
`config.yaml`'s antijam section was updated to the new schema (parses, and
`null_grid_deg`/`null_rank_cap` set to plausible values) but not re-tuned.
Re-running and re-tuning that campaign is the natural next step before
trusting the P1–P6 gate numbers again.

### P8 — Mode C: MUSIC DoA + predictive (anticipatory) nulling (amendment) — Status: **in-progress** (2026-07-27: on/off-first path IMPLEMENTED + demo; gates pass — see below. Drift/CV-Kalman follow-up not started)
Triggered by the customer's review request: exploit the jammer's **temporal
predictability** (on/off duty cycle; later, linear trajectory) so the null is
pre-positioned *before* the jammer returns/arrives, instead of reacting with the
`~1/(1−λ)` covariance lag. Targets the residual gaps on the non-stationary
scenarios (S5 toggling, S6 lifecycle; S2/S3 drift are a follow-up — see scope
note). Runs as a **new Mode C algorithm alongside** the P2 `lcmv` tracker, which
stays untouched as the reactive baseline for the comparison.

**In-session design decisions (2026-07-27):**
1. *Exploit which pattern first:* **on/off periodicity** (S5/S6). CV-Kalman drift
   prediction (S2/S3) is a follow-up on the same substrate.
2. *Predicted state → beamformer:* **hard LCMV null constraint at θ̂_pred**
   (`adapt_lcmv_null`), not a synthesized predicted-R MVDR. The constraint is
   deterministic in `w`, so it pre-nulls a currently-silent angle without needing
   energy there — exactly the pre-form-before-turn-on behavior.
3. *MUSIC model order:* **fixed at 2** (desired + 1 jammer, per locked scope).
   Jammer presence from the signal/noise **eigengap**, not order estimation.
4. *Coexistence:* **add** as algorithm `'predict'`; keep P2 `'lcmv'`.

**New modules** (flat `antijam_utils/`, R2020a base-MATLAB only — `eig`, `fft`; no
Signal Processing / Statistics Toolbox):
- `adapt_music_doa.m` — pure estimator: `R̂ → eig` → jammer `(θ̂,φ̂)`, presence
  (eigengap), power estimate, and (optionally) the pseudospectrum for the waterfall.
- `adapt_predict_init.m` / `adapt_predict_update.m` — reuses the P2 forgetting
  buffer for `R̂`; ring-buffers the presence indicator; `fft` periodogram detects
  the duty period + phase; schedules a pre-null `lead_steps` ahead of the predicted
  turn-on. Common `_init`/`_update` signature, consumes `obs.snapshots` ONLY.
- `adapt_lcmv_null.m` — multi-constraint LCMV, `C=[e_s, e_null]`, `g=[1;0]`
  (generalizes `adapt_lcmv`, which stays as-is). `e_null=[]` → plain MVDR.
- `plot_doa_waterfall.m` — the showcase 3-panel figure (below).

**Small edits to existing new-milestone files** (not `matlab_utils/`): add a
`'predict'` case to `closed_loop_run` (Mode C) and record DoA/presence/pspec
diagnostics into its `log`; add `'predict'` to `antijam.algorithms` at
implementation time (kept out of the list while only stubs exist so the campaign
doesn't error).

**Config schema addition** (`adapt.predict` block, parse-verified in this phase):
`presence_gap_db`, `buffer_len`, `min_periods`, `lead_steps`, `doa_stride`.
Model order is hardcoded to 2 (not a key). Values are initial guesses; expect
load-bearing keys to shift under P8 tuning (as with `adapt.spsa` in P3).

**Demo / showcase** (`plot_doa_waterfall`, reusing the P6 run/KPI harness):
(1) **DoA waterfall** — MUSIC pseudospectrum vs time, true `θ_j(t)` overlaid,
OFF gaps visible; (2) **presence periodogram** with the detected duty period
marked; (3) **SINR timeline** overlaying `lcmv` vs `predict` vs `oracle`, showing
`predict` pre-filling the recovery notch at each toggle.

**Success criteria (to gate at implementation):**
- On the toggling scenario (S5) and lifecycle (S6): recovery time after a
  turn-on ≤ the P2 `lcmv` baseline, with a clear reduction once the period is
  learned (target: near-zero re-acquisition lag after `min_periods` cycles).
- SINR availability on S5/S6 ≥ the P2 `lcmv` baseline (no regression); MUSIC
  DoA-RMSE reported (new KPI, MUSIC-derived) while the jammer is ON.
- Periodogram recovers the configured `toggle_period_s` within one grid bin on S5.
- `adapt_lcmv_null` reproduces `adapt_lcmv` exactly when `e_null=[]` (unit test).

**Scope note (honest, for the review):** on/off-first fully addresses S6 (static
angle held through the OFF gap — note S6 is a single burst, so the *periodogram*
finds no line; the "persist last predicted null" fallback carries it, not the
frequency detector). **S5 improves but will not fully close**: the angle drifts
*during* the OFF gap, so a held-static null is stale by `drift-rate × off-duration`.
Fully closing S5 needs the CV-Kalman drift predictor on this same substrate —
the deliberate P8 follow-up.

**New KPIs** (added to `kpi_evaluate`, existing `matlab_utils` metrics untouched):
DoA-RMSE (`θ̂_j` vs truth while ON) and prediction lead-time benefit on recovery.

**Implemented (2026-07-27):** all five files real + `closed_loop_run` `'predict'`
case + `test_antijam_predict.m` (4 gates pass: `adapt_lcmv_null`≡`adapt_lcmv`
when `e_null=[]` to 1e-12 and −60 dB constrained null; MUSIC angle+presence on a
toy analytic-R; on/off no-regression + period recovery + DoA-RMSE; Mode-C
contract). **Key implementation findings:** (1) MUSIC must use the
**normalized** metric `‖a‖²/‖Eₙᴴa‖²` — measured element-pattern column norms
vary strongly across the grid, so the unnormalized `1/‖Eₙᴴa‖²` spuriously peaks
at low-norm endfire directions (a ULA's constant `√N` norm hides this — the toy
tests passed while the real array failed). (2) During ON, let the measured `R̂`
form the null (identical to `lcmv`); use the explicit hard null ONLY for
pre-nulling — this makes `predict` never worse than `lcmv` and avoids
mis-steering on a transient estimate. (3) `buffer_len` must exceed
`min_periods × period_in_steps` or the periodogram never trusts a line; sub-bin
parabolic interpolation of the FFT peak is needed to time the pre-null (raw bin
was off by ~5 steps at period 200). (4) `adapt_lcmv_null` uses `pinv` on the
constraint gram to degrade gracefully if a predicted null coincides with the
steer direction. New required key `adapt.predict` (buffer widened to 1024,
`lead_steps` 6). **Demo** `scripts/run_mode_c_demo_script.m` on ManyDipoles/Theta,
on/off jammer (5 s off-gap, 10 s period, 100 s): oracle/lcmv/predict — dead time
**0.30 s (predict) vs 0.70 s (lcmv)**, availability 99.7 vs 99.3%, oracle gap
0.60 vs 0.78 dB, recovery 0.30 vs 1.22 steps; MUSIC DoA-RMSE ≈0° while ON;
periodogram recovers the 10 s period. Artifacts: `mode_c_comparison_*.png`,
`doa_diag_*.png`, `vid_*_{oracle,lcmv,predict}.mp4`, `mode_c_stats.txt`.

**Campaign integration (2026-07-27):** `'predict'` added to
`antijam.algorithms` in `config.yaml` — `run_antijam` is algorithm-agnostic
(loops `algorithms` → `closed_loop_run`), so no driver change was needed;
`kpi_evaluate` / `write_kpi_table` / `plot_antijam_report` all handle a
`predict` run (the extra DoA-diagnostic log fields are ignored by the generic
report paths). Smoke campaign (oracle/lcmv/predict × S1/S5/S6, 1 seed):
S1 (always-on static) predict ≡ lcmv exactly (no-regression invariant holds);
S5 (drift+on/off) predict avail 99.5 vs 99.2%, recovery 0 vs 1, oracle gap
0.86 vs 1.02, gain penalty −0.20 vs −0.32 dB; S6 (lifecycle) oracle gap 0.55
vs 0.75, gain penalty −0.10 vs −0.30 dB (predict restores the quiescent beam
better after turn-off). `plot_antijam_report` also given a light theme (white
figure/axes/legend via scoped `groot` defaults + `exportgraphics` white
background) so its figures match the demo — the report `print` path was
capturing the session's dark theme.

**Total-pol MVDR fix (2026-07-27):** `adapt_lcmv`/`adapt_lcmv_null` previously
passed the 2-column total-pol steering `e_s = [copol, cross]` straight into a
distortionless-on-EVERY-component LCMV (`wᴴe_copol = wᴴe_cross = 1`). That is
over-constrained for a rank-2 desired signal (independent per-component
amplitudes) and self-limiting when the components differ in magnitude — on
spacing0.6 at (θ30,φ120) the "oracle" hit only 14.3 dB SINR while 32.8 dB was
achievable (cross-pol there is ~18 dB stronger than co-pol, so forcing both to
unit gain wrecked the beam). Fixed with new `adapt_maxsinr.m` (rank-r max-SINR =
principal generalized eigenvector of `(R_s, R+load·I)`, optional exact hard null
via projection onto `{w : e_nullᴴw = 0}`); `adapt_lcmv`/`adapt_lcmv_null`
delegate to it for `n_c ≥ 2` and keep the exact rank-1 MVDR/LCMV for single
components (byte-identical, `max|Δw| = 0`; all suites still pass). Verified: the
`total` oracle now equals the true SINR upper bound (32.8 dB jammer-off,
30.2 dB with a J/N 20 dB jammer at (30,150), −77.6 dB null; target directivity
0.3→15.8 dBi). New regression gate `test_antijam_predict` gate 5. This fixes the
oracle/`lcmv`/`predict` for `total` pol; single-component pol was never affected.

**Not validated yet:** the **drift / CV-Kalman** predictor (S2/S3) — the P8
follow-up. `plot_doa_waterfall` implemented as a DoA-track + periodogram figure
(not a literal pseudospectrum waterfall — 2-D geometry). Also inherits the P7.4
caveat: the full 2-D P1–P6 Monte-Carlo campaign has not been re-run/re-tuned, so
those KPI numbers remain stale (the smoke run above is 1 seed, 3 scenarios).

**Covariance-tracker batching fix + re-tune (2026-08-01):** `adapt_tracking_update`/
`adapt_predict_update` were folding the `K` snapshots/step in via K sequential
per-column `(1-lambda)` recursions instead of one batch-averaged per-step
update — an accidental `lambda^K` per-step decay that the original P2 sweep's
`forgetting_lambda = 0.98` was unknowingly tuned against. Fixed to match the
plan's own formula (Section 2); re-swept `forgetting_lambda` to **0.90** on
both `test_antijam_tracking` and `test_antijam_predict` (full margin 0.85–0.92).
`mode_c_demo` steady-state oracle gap improved 8.64→5.65 dB (`lcmv`), 5.78→4.74
dB (`predict`); residual gap on that specific demo scenario is attributed to a
genuinely sharp optimum (embedded cross-pol gain at the demo's jammer angle
rivals the desired-signal eigenvalue), not the tracker itself — see
`docs/notes.md` `[P8]` 2026-08-01 entry for full detail. This changes the
effective memory time constant of every P2/P8-gated result quoted above in this
section (all measured under the old `lambda=0.98`/buggy-decay combination);
none of those suite gates were re-run against the fix beyond the two files
named, so treat quantitative KPI numbers elsewhere in this P8 section (and P2's
1-D-cut-era numbers) as superseded pending a full re-run.

**Weight-vector smoothing feature (2026-08-01/02):** even with the batching
fix, `lcmv`/`predict` still "snap" to a fresh closed-form solution every step —
no amount of `R_hat` smoothing produces a genuine gradual approach to the
optimum on a sharp-optimum scenario. Added optional `adapt.weight_smoothing_mu`
(new key, default off): phase-aligned EMA on the *applied* weights,
`w <- (1-mu)*w + mu*w_target`, in `adapt_tracking_update.m`/
`adapt_predict_update.m`. Opt-in only (`mu` absent/1 reproduces prior behavior
exactly, verified `max|Δw|=0`); costs measurable availability/reacquisition
speed as mu drops, so it has NOT been folded into the P2/P8 gate suites (which
still test the un-smoothed tracker) — only `config.yaml`'s demo/campaign
tuning uses it (`weight_smoothing_mu: 0.25`, alongside `snapshots_per_step:
128`). Final mode_c_demo (25 dB jammer, 20 s toggle scenario): avail
96.1%/96.1%, oracle gap 3.06/2.93 dB (`lcmv`/`predict`), visibly smooth climb
after every turn-on. See `docs/notes.md` `[P8]` 2026-08-01/02 entries for the
full sweep data and the base-tuning correction (reverted `diagonal_loading_db`/
`forgetting_lambda` to 10 dB/0.90 after re-sweeping against the demo's actual
scenario — an earlier "14 dB/0.95" combo had been fit to a since-changed,
much weaker 1 dB jammer scenario).

---

## 5. KPIs and standard scenario suite

### KPIs (implemented in `kpi_evaluate.m` — existing `matlab_utils` metrics untouched)
1. **SINR availability** *(headline)* — fraction of timeline with SINR ≥ `sinr_min_db`.
2. **Recovery time** — updates/probes from a jammer event (jump, power step, turn-on) until SINR re-crosses threshold.
3. **Null-pointing error** — |achieved null angle − θ_j(t)|, tracked over time.
4. **Peak-gain penalty** — main-beam gain loss vs jammer-free Milestone-1 optimum (dB).
5. **Oracle gap** — SINR shortfall vs perfect-knowledge LCMV (dB), per instant and averaged.

### Standard scenario suite (defined in config; seeds fixed)
| ID | θ_j behavior | Power | Notes |
|---|---|---|---|
| S1 | Static, random angle per seed | constant, J/N = 20 dB | identification baseline |
| S2 | Linear drift 2°/s | constant | primary tracking case |
| S3 | Drift + 10° jump mid-run | constant | recovery-time case |
| S4 | Static | power step +10 dB mid-run | robustness |
| S5 | Drift, on/off toggling (duty 50%) | constant | non-stationarity stress |
| S6 | Static; silent 0–60 s, ON 60–120 s, silent 120–180 s | `window` mode | null lifecycle: null forms on turn-on, quiescent beam restored after turn-off (added 2026-07-19 per Snir; `tests/test_antijam_lifecycle.m`) |

Angles avoiding the main-beam sector by a configurable guard (jammer inside the main beam is out of scope for a single-jammer milestone — note as known limitation).

---

## 6. Config schema additions (sketch)

> **P0 note:** the schema below is now implemented in `config.yaml` (authoritative copy — it additionally carries `antijam.cut_type` / `cut_theta_deg` / `cut_phi_deg`, `antijam.algorithms`, and the full `antijam.scenarios` suite). Parse verified with `read_config_yaml` on 2026-07-18.

Parsed by the existing `read_config_yaml` — which supports nested **block** mappings, inline flow *lists* `[a, b]`, and block sequences, but **not** inline flow mappings (`{a: 1, b: 2}`). All nesting below is therefore block-style.

```yaml
antijam:
  theta_s_deg: 0.0
  sinr_min_db: 10.0
  sigma_s_db: 0.0          # signal power ref
  jn_ratio_db: 20.0        # jammer-to-noise
  guard_deg: 10.0          # jammer exclusion around main beam
sim:
  dt_s: 0.05               # closed-loop step
  duration_s: 60.0
  snapshots_per_step: 16   # Mode C
  seed: 1234
adapt:
  diagonal_loading_db: -10  # relative to noise floor
  forgetting_lambda: 0.98
  spsa:            # a_k = a/(A+k)^alpha, c_k = c/k^gamma, step norm <= step_max
    a: 2.0
    c: 0.2
    alpha: 0.602
    gamma: 0.101
    A: 15
    step_max: 0.3
agent:
  null_grid_deg: 5.0        # codebook arm spacing
  null_width_deg: 10.0
  method: thompson          # thompson | swucb
  discount: 0.97
  window: 50                # for swucb
```

---

## 7. Open decisions — RESOLVED in P0 (session 2026-07-18)

1. **`sinr_min_db` value** — ✔ default 10 dB adopted (customer may revise; plain config key).
2. **Codebook arm spacing** — ✔ default 5° adopted (`agent.null_grid_deg`).
3. **Reward definition for bandit** — ✔ SINR in dB, clipped at [−10, +40].
4. **Field convention for the SINR engine** — ✔ **config-selectable**, mirroring `run_optimization`'s `polarization` handling (`copol` / any named component / `total`). The convention enters the engine only through whether a secondary component cut is supplied (`E_secondary` empty vs not); `total` uses incoherent power sum with equal power split between components (unpolarized-source assumption — verify in P1 tests). Codebook arms are synthesized under the same `polarization` setting.
5. **Snapshot count per step (Mode C)** — ✔ default 16 (`sim.snapshots_per_step`).

---

## 8. MATLAB implementation constraints & validation

**Compatibility target: R2020a + Optimization Toolbox only** (same constraint as the Milestone-1 compatibility pass — see `docs/MATLAB_R2020a_changes.md`). No Phased Array, Signal Processing, or Statistics Toolbox. No post-R2020a syntax/functions.

- **Repo conventions apply:** flat file = function in `antijam_utils/`, module prefix in the filename, header-comment docstrings in the existing port style, 1-based indexing throughout.
- **No classes** in `adapt_*` / `agent_*` / `sim_*` — plain functions + state structs, `<name>_init` / `<name>_update` pairs.
- **Determinism:** every stochastic component seeds via `rng(seed, 'twister')` from config; all P1/P5 criteria measured under fixed seeds.
- **Complex Gaussian sampling:** `(randn(...) + 1i*randn(...))/sqrt(2)` (base MATLAB; no toolbox needed).
- **Thompson sampling:** rewards are real scalars — `randn`-based Gaussian posterior sampling only.
- The Python-era `# [MATLAB]` porting-flag convention does **not** apply here — this code is MATLAB-native.

**Validation strategy (no Python reference exists for this milestone):** unit tests in `MATLAB/tests/test_antijam_*.m` against **analytic and toy-case ground truth**, in the style of the existing test suite but with fixtures computed in-test rather than loaded from Python-generated JSON:
- closed-form oracle LCMV null depth and SINR on small arrays,
- hand-computed 2-element SINR (1e-10 relative),
- covariance-convergence rates vs theory,
- fixed-seed regression on scenario trajectories,
- codebook arm gate checks (P4 criteria) as automated assertions.

---

## 9. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| Bandit tuning (discount vs drift rate) takes longer than planned | Medium | P3 SPSA baseline is a shippable fallback; sweep discount in P5 with automated grid |
| Codebook holes between arms (jammer between sectors) | Low–Med | Coverage plot in P4 gate; overlap arms via `null_width_deg` > spacing |
| Covariance mode over-nulls with limited snapshots | Medium | Diagonal loading swept in P2; loading level in config |
| Pattern-level sim judged insufficient by customer | Low | Interfaces designed so a signal-level engine can replace `sim_engine_step` without touching algorithms |
| `read_config_yaml` subset can't express a needed config shape | Low | Schema in Section 6 already restricted to the supported subset; verified in P0 before implementation |
| R2020a constraint blocks a convenient built-in | Low | Reimplement from definition, as done for the window functions in the Milestone-1 port |

---

## 10. Session log pointer

Per project convention, every working session on this milestone appends an entry to `docs/notes.md` referencing the phase ID (e.g., `[P2] tuned lambda sweep, chose 0.98 — oracle gap 0.7 dB on S2`). Phase `Status` fields in Section 4 are the single source of truth for progress.
