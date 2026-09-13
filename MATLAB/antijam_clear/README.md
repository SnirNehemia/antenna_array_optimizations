# `antijam_clear` — a readable anti-jam beamformer

A deliberately simple reimplementation of the anti-jam milestone, built to be
**explained** rather than to win. Two approaches to the same problem, kept
separate on purpose:

1. **Detect, then null** — estimate where the jammer is and what it is doing,
   then shape the pattern accordingly. Explicit, geometric, legible.
2. **Maximise SINR directly** — solve in closed form for the weights that
   maximise output SINR, never locating the jammer at all.

It depends on `MATLAB/matlab_utils/` for CST parsing and pattern maths, and on
nothing in `MATLAB/antijam_utils/`. Results are written to
`MATLAB/antijam_clear/results/`, separate from the repository-level `results/`.

Runs on MATLAB R2020a. Background and the settled facts this is built on:
`docs/antijam_rewrite_brief.md`.

---

## Deliberate omissions

Everything here was present in the previous implementation and is left out on
purpose. Each entry says what it cost, so the trade is on the record and so the
delivery version can pick any of them back up without re-deriving why.

### Fast presence covariance — *flagged for the delivery version*

A **second** covariance EMA at a much shorter memory (`lambda = 0.5`), used
**only** for presence detection and never to form weights. The principle: a
detector slower than the beamformer cannot resolve edges the beamformer already
smooths.

- **Measured worth:** presence error against the true duty cycle fell from
  0.46 / 0.32 / 0.14 to 0.13 / 0.06 / 0.03 at toggle periods of 4 / 10 / 25 s.
  The anticipatory branch went from firing on ≤1.3% of steps to 1.4–4.7%.
- **Cost of omitting it — measured here, and worse than "one horizon late".**
  The source threshold sits just above the *noise floor*, but a strong jammer's
  eigenvalue starts far above it: on `spacing0.6` at JNR 20 dB it peaks **23.7 dB**
  above the threshold. The EMA decays it by only `−10·log10(λ) = 0.46 dB/step`,
  so the residual needs **≈52 steps** to fall below the threshold. The jammer is
  off for 10 steps at a time, so **the off state is never observed at all** —
  presence reads 100% against a true duty cycle of 50%, at every JNR from 3 to
  20 dB. The general form is

  `steps to notice the jammer left = (jammer eigenvalue excess in dB) / (0.46 dB per step)`

  which is set by *jammer strength*, not by the toggle period.
- **What it does not cost:** nothing, in this implementation. `steady` and
  `onoff` both hold the null at the last known position, so the *action* is
  identical and only the reported label is wrong. The fast covariance is needed
  for **anticipation** — acting before the jammer returns — not for basic
  nulling. That is the capability being deferred.
- **How to add it back:** call `sample_covariance` a second time with its own
  `previous_covariance` and its own lambda, and feed it only to the presence
  test in `classify_jammer_motion`. One extra line of carried state; it must
  never reach the weight solve.

### CV-Kalman direction predictor

A constant-velocity Kalman filter on `[θ, θ̇, φ, φ̇]` extrapolating a derived
lead of 14 steps (10 covariance horizon + 4 application and smoothing lag).

- **Measured worth:** +35 points on the drifting-jammer score (55.0 → 90.1),
  and provably inert on steady and on/off runs.
- **What replaces it here:** a one-line arithmetic lead. The beamscan lags the
  jammer by `rate × (mean age of the covariance)`, so the null is aimed that far
  ahead, using the fitted rate and a quantity that follows from λ. Nothing is
  tuned to performance. On the demo scenario it takes the aim error from −4.50°
  to −0.37°, mean |error| from 4.65° to 0.75°, and the score from 65 to 99.
- **Which horizon — worth about 1 dB.** The covariance weights a block k steps
  old by `(1−λ)λᵏ`, and that weighting has two one-number summaries:
  `1/(1−λ) = 10` (effective window length, sets memory and the motion window)
  and `λ/(1−λ) = 9` (mean age of the data, sets the lag). They differ by exactly
  one step. **The lead must use the mean age.** Measured: dividing the observed
  lag by the drift rate gives 9.00 steps exactly at 0.10 and 0.20°/step. Using
  the window length instead over-leads by one step's motion every step and costs
  0.4–1.3 dB across the qualified range.
- **Cost of omitting the real filter — the drift envelope.** The lead multiplies
  the fitted rate by 9, so it amplifies any error in that rate ninefold. It pays
  while the rate is well measured and stops paying when it is not:

  | drift rate | detect+null | achievable | aim error | score |
  |---|---|---|---|---|
  | 0.25°/step | 11.1 dB | 11.2 dB | 0.34° | 98 |
  | 0.50°/step | 8.9 dB | 11.2 dB | 0.72° | 70 |
  | 1.00°/step | 1.9 dB | 11.1 dB | 3.75° | 25 |
  | 1.50°/step | −4.0 dB | 11.0 dB | 5.61° | 3 |

  (90° traverse; `results/drift_envelope.csv`.) **This version is qualified to
  roughly 0.5°/s and broken by 1.5°/s.** The previous implementation, with the
  CV-Kalman filter, was qualified to about 10°/s — an envelope some **20×
  wider**. That is the honest price of leaving the filter out, and it is a
  narrower *envelope*, not a few dB.
- **The lead is capped at the guard sector** (half a beamwidth), because beyond
  that there is no measurement bearing on where the jammer will be. The cap
  reuses `guard_deg` and adds no constant. It is inert below 1.26°/step and
  rescues the high-rate collapse (2°/step: −14.1 → −4.7 dB).
- **Measured: adding the filter back is NOT the way to widen the envelope.**
  The aim is `angle + rate × mean age`, and a filter can only improve the *rate*.
  Feeding the loop the **true** rate bounds any rate estimator that could exist:

  | drift rate | ships | perfect rate | perfect angle | achievable |
  |---|---|---|---|---|
  | 0.25°/step | 11.1 dB | 11.1 dB | 11.2 dB | 11.2 dB |
  | 0.50°/step | 8.9 dB | 10.4 dB | 11.2 dB | 11.2 dB |
  | 1.00°/step | 1.9 dB | 3.1 dB | **11.1 dB** | 11.1 dB |
  | 1.50°/step | −4.0 dB | −4.2 dB | **11.0 dB** | 11.0 dB |

  (`results/predictor_bound.csv`.) **A perfect rate estimator is worth ≤1.5 dB**
  — and above 1°/step it is *worse*, because past a beamwidth of travel per
  horizon the lag model stops holding and the correction overshoots.
  **Eight dB sits in the angle, at most 1.5 dB in the rate.** Given the right
  angle the two-point null already matches the oracle to 0.1–0.3 dB.
- **The real lever is λ.** Shortening the memory is worth **+8.3 dB at 1°/step**
  (λ 0.90 → 0.70; 5 seeds, spread ±0.1–0.4 dB) — but it **breaks the on/off case**
  (12.1 → 4.4 dB): a short memory loses the jammer while it is off, the beamscan
  peak wanders onto noise (out to 143°), presence detection is too late to catch
  it, and the motion classifier then reads that noise as drift and amplifies it
  by the lead. **The fast presence covariance above is exactly what makes a short
  memory safe.** Order of work: presence detector first, shorter memory second.
  The CV-Kalman filter is not on the path.

### Two-component ("total") polarization

Only one polarization component is loaded. With two, each source becomes rank-2
instead of rank-1, MUSIC needs a model order of two per source, and every
`|wᴴe|²` becomes a sum over components.

- **Cost of omitting it:** results are co-pol (or θ-pol) only, and do not
  account for a jammer arriving on the orthogonal polarization.
- **How to add it back:** call `load_array` twice and carry a second cube; see
  the `[FUTURE]` notes in `load_array`, `steering_vector` and
  `simulate_snapshots` for the exact shape changes.

### Mode S, the frequency/notch layer, campaign runners

Scalar-SINR algorithms (SPSA, bandit), the waveform and RF-notch layer, and
every campaign runner, arm sweep and opt-in config block. Out of scope for a
readable Mode C implementation; none of it is load-bearing for either approach
here.

---

## Both direction finders are kept, and the simpler one is used

`estimate_jammer_angle` computes two spectra every step and reports both. The
**beamscan** places the null; **MUSIC** is reported for comparison and never
used. There is deliberately no switch between them — the comparison is the
output, not a setting.

MUSIC is the textbook choice and is sharper on a stationary source, because it
measures a null in the noise subspace rather than a peak. On a **moving** source
it is measurably worse, and the reason is worth knowing: a jammer smeared across
ten steps of covariance memory is not a point source — the source count reads 3
on 91 of 120 steps of a drifting run — so MUSIC makes its noise subspace
orthogonal to the jammer's whole *track*, its spectrum goes near-singular along
a broad arc, and the peak inside that arc is numerically arbitrary.

| estimator | median lag | mean \|lag\| | max \|lag\| | behaviour |
|---|---|---|---|---|
| MUSIC | 3.50° | 5.03° | **13.50°** | sticks, then jumps |
| beamscan | 4.50° | 4.65° | **6.00°** | smooth |

The beamscan's lag is exactly `rate × mean age`, because a power measure peaks at
the power-weighted centre of where the jammer has been. That is a *predictable*
error, and the lead correction removes it; MUSIC's 13.5° excursion cannot be led
out. Switching to the beamscan took the drift case from 7.5 dB / score 66 to
11.1 dB / score 99 against an achievable 11.5, measured with identical settings
and identical noise, changing only which spectrum feeds the null.

On a stationary jammer the two agree to 0.00°. MUSIC's real advantage is
resolving two closely-spaced sources, which a single-jammer problem never asks
for — and the guard sector already excludes the region near the signal.

---

## Compute cost — the scan dominates, and Approach 2 avoids it entirely

Measured per adaptation step (MATLAB R2020a, one core, double precision,
16-element array, 181×360 grid). Relative costs are structural; absolute numbers
would improve a lot in an embedded single-precision implementation.

| per step | ms | share |
|---|---|---|
| `sample_covariance` | 0.055 | 0.2% |
| `estimate_jammer_angle` (both scans) | **30.415** | **95.9%** |
| — beamscan over the grid | 15.001 | 47.3% |
| — MUSIC over the grid | 10.883 | 34.3% |
| — re-normalising the manifold | 5.213 | 16.4% |
| `detect_and_null` (2×2 solve) | 0.065 | 0.2% |
| `max_sinr_weights` (no scan) | 0.057 | 0.2% |
| **full step as shipped** | **31.7** | 100% |

**The search is ~4,000× the arithmetic of the solve**: a full 2-D scan is about
133 MFLOP and streams 16.7 MB of manifold, against ~33 kFLOP for the entire
max-SINR solution. Direction finding is 96% of the per-step cost.

| variant | ms/step | max rate |
|---|---|---|
| as shipped | 31.7 | 32 Hz |
| beamscan only, manifold precomputed | 14.7 | 68 Hz |
| coarse 5° grid then ±5° refine at 1° | **0.46** | **2,170 Hz** |
| 1-D θ scan at a known azimuth | 0.075 | 13,300 Hz |
| **Approach 2, no scan at all** | **0.062** | **16,200 Hz** |

Three reductions are free and not yet applied:

- **Precompute the unit manifold** in `make_array` — it is a fixed property of
  the array, currently rebuilt every step. 16% for nothing.
- **Drop MUSIC from the hot path** — it is reported-only and never places a
  null. Another 34%.
- **Coarse-to-fine** — 69× cheaper and *exact* here: identical peak on all 180
  steps of the steady and drift runs, 0° difference.

**Structural point:** Approach 2 needs no scan at all — no grid, no search, no
resolution/compute trade. If compute is the binding constraint it is the
deployable algorithm, and its 3.4 dB deficit on drift must be weighed against
being 500× cheaper than the method that beats it.

---

## Not to be reintroduced

Two features were built, measured at scale and rejected. They are recorded here
so they are not proposed again as improvements.

- **Graded / confidence-weighted null release.** Scored 76.1 against the binary
  policy's 80.6 and was worse on half the cells. More decisively, *best-of-
  {none, binary} per cell* — a bound no real switching rule can beat — scores
  81.5. The entire family is capped at **+0.9 mean and zero extra passing
  cells.**
- **Data-driven adaptive diagonal loading.** Measured at **+0.00 pp across all
  90 cells**, and it was the direct cause of two arrays being unable to run.
  (Diagonal loading itself is mandatory — it is the *adaptive* variant that is
  worthless.)
