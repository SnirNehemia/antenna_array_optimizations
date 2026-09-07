# Mode C test campaign — running findings

## F1. ManyDipoles has an EXACT theta-mirror ambiguity (new, 2026-09-07)
`e(theta,phi) == e(180-theta,phi)` to 2.8e-6 relative over the whole stack.
Steering coherence = 1.0000 exactly at every probed angle.

- This is NOT the antipodal pair (180-theta, phi+180) whose low coherence
  (0.013-0.30) the 2026-09-06 note used to RETRACT the front/back diagnosis.
  That retraction measured the wrong pair; the ambiguity is real, it is just
  a different one.
- Physically expected: a planar array in the z=0 plane with a common phase
  centre and no geometric phase re-added is exactly up/down degenerate.
  The other two arrays have ground planes (patches/monopoles) which break it:
  mirror coherence 0.85 -> 0.08 (patchs) and 0.87 -> 0.53 (spacing0.6) as
  theta goes 110 -> 170.
- CONSEQUENCE A (benign): because the two steering vectors are IDENTICAL, a
  null placed at the mirror angle IS the null at the true angle. `lcmv` never
  uses an angle, and `predict`'s hard null is therefore also unharmed. The
  60 deg "DoA error" measured below is a LABELLING artifact, not a
  performance defect.
- CONSEQUENCE B (harmful): any DoA-error KPI on this array is meaningless
  unless computed modulo the fold, and a naive CV-Kalman on theta WILL break
  -- the MUSIC peak flips between branches, producing a garbage velocity.
  The CV predictor must be fold-invariant.
- CONSEQUENCE C (deliverable): the array as exported cannot tell up from
  down. Worth stating to the customer.

Measured DoA error, DRIFT (2 deg/s in theta), 30 s, stride 1, presence 100%:
| array | median | p90 | max | RMSE | >20 deg |
|---|---|---|---|---|---|
| patchs_with_monopoles | 1.00 | 1.40 | 1.50 | 1.03 | 0.0% |
| spacing0.6 | 1.00 | 1.40 | 1.50 | 1.06 | 0.0% |
| ManyDipoles | 60.20 | 105.20 | 120.00 | 67.96 | 83.4% |

## F2. `predict` cost, measured (10 s DRIFT run, x lcmv)
| array | stride 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| patchs_with_monopoles | 42.5x | 25.8x | 16.4x | 11.3x |
| spacing0.6 | 126.3x | 71.5x | 41.8x | 28.0x |
| ManyDipoles | 4.6x | 3.5x | 2.8x | 2.5x |

DoA RMSE cost of striding (1-deg arrays): 1.00 -> 1.45 (stride 4) ->
2.30 deg (stride 8) on spacing0.6. Campaign runs at stride 1.

## F0. Gate baseline (2026-09-07)
79/80 tests pass. Sole failure `test_metrics/test_evaluate_metrics_matches_python`
is the known pre-existing MATLAB-vs-Python phi-wrap discrepancy in
matlab_utils/, unrelated to anti-jam. All 23 anti-jam gates green.

## F3. The DRIFT failure IS angle lag, and it is exactly one lambda-horizon
Decisive experiment (patchs_with_monopoles, sep45_th, DRIFT 2 deg/s, 30 s,
seed 1234; hard null steered from scenario TRUTH -- diagnostic only):

| variant | track% | mean gap dB |
|---|---|---|
| oracle | 100.0 | 0.00 |
| lcmv (shipped) | 44.8 | 4.00 |
| clairI (null @ true CURRENT angle, R = eye) | 94.2 | 0.89 |
| clairR (null @ true CURRENT angle, R = R_hat) | 94.8 | 0.88 |
| null @ truth delayed 5 steps | 73.5 | 2.28 |
| null @ truth delayed 10 steps | 49.8 | 3.88 |
| null @ truth delayed 15 steps | 38.6 | 5.29 |
| null @ truth delayed 20 steps | 27.3 | 6.41 |

Readings:
1. lcmv is reproduced by a truth-steered null delayed 10-11 steps
   (49.8% / 3.88 dB vs lcmv's 44.8% / 4.00 dB). 10 steps = 1/(1-lambda) at
   lambda = 0.90 -- the covariance forgetting horizon, exactly. At dt = 0.05 s
   and 2 deg/s that is 0.5 s = 1.0 deg of angular lag.
2. Perfect CURRENT-angle knowledge lifts the score 44.8 -> 94.8%. So the drift
   gap is an ANGLE-LAG problem and a predictor can close it. CV-Kalman is
   justified on measurement, not on assumption.
3. clairI ~= clairR: once the angle is known the covariance contributes
   essentially nothing. The predictor's job is angle, not covariance.
4. CEILING is ~95%, not 100% -- residual is grid quantization, finite null
   depth and the weight_smoothing_mu = 0.25 lag. A realistic CV target is
   90-95%, i.e. just over the pass mark, not a blowout.

## F4. CV-Kalman predictor: implemented, and the lead is DERIVED
New `adapt_cv_init` / `adapt_cv_update`, wired into `adapt_predict_update` as
an opt-in `adapt.predict.cv` block. With the block absent the code path is
byte-identical to pre-P12b.

**The lag budget, measured.** The raw MUSIC angle best matches truth delayed
exactly 10 steps (mean |err| 0.27 deg at tau=10, rising either side) -- because
MUSIC is computed from the SAME exponentially-forgotten R_hat, so it inherits
the 1/(1-lambda) = 10-step horizon. Add the 1-step application delay and the
weight_smoothing_mu = 0.25 first-order lag (~4 steps) and the total is 14.
`lead_steps = 14` is therefore derived, and the design rule is that it should
track 1/(1-lambda) if lambda ever moves.

**Lead sweep, patchs_with_monopoles, track% (lcmv baseline in last column):**
| drift | L=10 | L=12 | L=14 | L=16 | L=18 | lcmv |
|---|---|---|---|---|---|---|
| 1 deg/s | 97.5 | 95.0 | 90.2 | 85.7 | 80.7 | 64.4 |
| 2 deg/s | 83.7 | 88.7 | 90.3 | 91.0 | 88.5 | 44.8 |
| 4 deg/s | 56.8 | 65.4 | 71.1 | 71.1 | 71.4 | 11.3 |
| MEAN | 79.3 | 83.0 | **83.9** | 82.6 | 80.2 | |

The optimum lead is rate-dependent (bias-variance: leading amplifies velocity
error). A single fixed value is a compromise -- the same structural honesty the
P12 lambda sweep reached. It is never WORSE than lcmv at any rate tested.

**Regression safety is exact, not incidental.** The min_speed_deg_s gate keeps
CV out of non-drifting runs, and STATIC / WINDOW are BIT-IDENTICAL with cv on
vs off (max |dSINR| < 1e-12).

**Drift + 10 deg jump (S3-like): 36.6% -> 86.7%.**

**Interaction found, NOT applied:** weight_smoothing_mu is itself a drift cost.
lcmv alone goes 44.8 -> 55.4 -> 59.1% as mu goes 0.25 -> 0.5 -> 1.0, and
CV+mu=1.0 reaches 93.0% at L=12. mu was tuned for other scenarios, so changing
a global default is out of scope here -- but it is a real second lever.

## F5. Calibration mismatch is the dominant deliverability risk (NEW capability)
The milestone plan listed steering-vector mismatch as "the largest untested
failure mode, and structurally untestable today" -- `sim_engine_step` computed
the signal power from the IDENTICAL `e_s` handed to the beamformer. P12b adds
an OPT-IN `assumed` argument to `closed_loop_run` (no `sim_` module touched):
the engine propagates the TRUE patterns, the algorithm is handed assumed ones.
Self-consistency verified: passing the true stacks as `assumed` reproduces the
no-assumed run to max |dSINR| = 0.00e+00.

**Robustness curve.** Per-element calibration error g_n = exp(j*e_n), e_n ~
N(0, sigma_phase). DRIFT 2 deg/s, patchs_with_monopoles, fixed 10 dB loading,
12 error draws per level, oracle-tracking score:

| phase err [deg RMS] | lcmv median (p10-p90) | predict+CV median (p10-p90) |
|---|---|---|
| 0.0 | 44.8 | 90.3 |
| 0.5 | 43.3 (36.8-46.3) | 89.1 (84.8-91.0) |
| 1.0 | 37.4 (19.4-44.7) | 85.4 (39.7-90.0) |
| 1.5 | 24.6 (4.0-40.8) | 72.9 (4.7-86.0) |
| 2.0 | 10.5 (0.5-35.0) | 45.2 (1.1-76.8) |
| 3.0 | 1.8 (0.3-15.4) | 6.8 (0.3-35.3) |
| 5.0 | 0.3 (0.3-1.6) | 0.6 (0.3-4.4) |

Readings:
1. **The usable region is below ~1 deg RMS phase error.** At 3 deg both
   algorithms are at the floor. This is the classic MPDR signal-cancellation
   sensitivity -- the snapshots contain the desired signal, so a mis-specified
   constraint direction makes the beamformer null its OWN signal.
2. **predict+CV keeps its advantage at every error level** (2-4x lcmv until
   both die). The worry that a HARD null at an ASSUMED steering column would be
   MORE fragile than a covariance-derived null is NOT supported.
3. **The spread is enormous.** At 2 deg the p10-p90 for predict+CV is 1.1-76.8:
   one array build could be fine and the next catastrophic. Any acceptance
   statement about calibration needs a distribution, not a point.
4. Adaptive (P9) vs fixed loading measured IDENTICAL under mismatch in this
   regime, consistent with config.yaml's note that loading_factor_db = 0
   reproduces the tuned ~10 dB here. (An earlier reading of mine that adaptive
   loading caused the collapse was wrong -- it was a different random error
   draw, not a different loading mode. Corrected here.)

**Consequence for the recommendation:** the CV predictor is a real ~2x gain on
drift, but it does NOT change the calibration exposure, which is the larger
risk to fielding. Calibration tolerance belongs in the spec.

**Not a calibration test:** data/spacing0.6_disturbed3 vs data/spacing0.6 is a
file-for-file pair but the two arrays differ grossly (steering coherence median
0.40, min 0.15), so every algorithm sits at 0.2%. It is a "wrong array" test,
not a "mis-calibrated array" test. The parametric per-element model above is
the right instrument.

## F6. The mirror fold auto-detects correctly and rescues the DoA (end-to-end)
`adapt_cv_init` measures each array's own mirror coherence at init:

| array | mirror coherence | fold enabled |
|---|---|---|
| ManyDipoles | 1.0000 | yes |
| patchs_with_monopoles | 0.0505 | no |

On ManyDipoles / DRIFT, the MUSIC DoA error against truth is:
- **raw: median 139.2 deg**, >20 deg on 95% of steps
- **modulo the fold: median 1.5 deg**, >20 deg on 18% of steps

So MUSIC was accurate the whole time -- it was reporting the mirror branch, and
the array's 5 deg grid explains the residual 1.5 deg. This closes F1: the "68
deg DoA RMSE" was entirely a labelling artifact, and the plan's open P8 item
"MUSIC DoA front-end untested on real element patterns" now has an answer --
it works, once the ambiguity is accounted for.

## F7. The one residual failing cell is GRID-limited, not algorithm-limited
`spacing0.6 / sep45_th / DRIFT` (theta_j 135 -> 165 deg, 16 elements) is the
only profile-A DRIFT cell still under 90% with CV on (84.4% in the campaign;
83.5% reproduced single-seed).

Isolated:
- DoA quality is NOT the cause: median error 1.00 deg, presence 100%.
- Lead mistuning is NOT the cause: the sweep is flat and still climbing at
  L = 24 (78.9 / 82.2 / 82.7 / 83.5 / 83.7 / 84.4 / 85.0 / 85.9 for
  L = 8/10/12/14/16/18/20/24) -- a further +2.4 pp for a lead 70% larger.
- The CEILING for this cell (truth-steered null, same loop) is **98.3%**.

So ~14 points sit between a perfectly-steered null and a null steered by a
1-deg-quantized DoA. This is the highest-DoF array in the suite (16 elements),
so its null is the narrowest, and residual pointing error costs the most there.
The binding limit is the 1 deg export grid, not the estimator.

**Concrete recommendation:** interpolate the MUSIC peak sub-grid (parabolic
interpolation on the pseudospectrum -- `adapt_predict_update` already does
exactly this for the FFT period estimate, so the technique is in-tree). Expected
to matter most on high-element-count arrays and to be nearly free.
