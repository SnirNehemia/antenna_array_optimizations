# Phase O — stationary on/off jammer, all arrays. Running findings.

## O-F1. Array inventory: 5 viable, 2 degenerate
| array | n_el | grid | pol | mirror coh | best dBi |
|---|---|---|---|---|---|
| patchs_with_monopoles | 6 | 1 deg | dual | 0.57 | 9.9 |
| spacing0.6 | 16 | 1 deg | dual | 0.53 | 19.0 |
| spacing0.6_disturbed3 | 16 | 1 deg | dual | 0.37 | 16.1 |
| Monopoles | 14 | 1 deg | dual | 0.32 | 16.9 |
| ManyDipoles | 20 | 5 deg | single | 1.0000 | 17.6 |
| patch_back2back | 2 | 1 deg | dual | 0.20 | 9.2 |
| Dipole | 1 | 5 deg | single | 1.0000 | 3.2 |

`spacing0.6_disturbed3` is a DIFFERENT array from `spacing0.6` (mutual steering
coherence median 0.40), not a mis-calibrated copy. Its potential is 2.9 dB below
`spacing0.6`'s, so it tests adaptation to a degraded aperture.

## O-F2. Target direction is itself a difficulty axis (never varied before)
Quiescent directivity toward candidate targets swings ~18 dB on one array
(spacing0.6: 0.72 dBi at (150,100) vs 18.63 dBi at (20,260)). HPBW ranges
24-120 deg across (array, target). Every campaign to date fixed the target at
(90, 260) -- which on ManyDipoles is the mirror-symmetry plane, i.e. that
array's single most favourable target.

**Consequence:** `guard_deg = 5` (derived once from ManyDipoles' 15 deg HPBW)
is meaningless against patchs' 120 deg beam at (90,260). The guard sector MUST
be derived per (array, target).

## O-F3. TWO defects stop the stack running at all on small arrays
Measured on `Dipole` (1 el) and `patch_back2back` (2 el, dual-pol):

| array | shipped config | with FIXED loading |
|---|---|---|
| Dipole / lcmv | ERROR adapt_tracking_init:TooFewElements | **100.0%** (SINR 5.92 = oracle 5.92) |
| Dipole / predict | ERROR adapt_predict_init:TooFewElements | ERROR adapt_music_doa:TooFewElements |
| patch_back2back / lcmv | ERROR adapt_tracking_init:TooFewElements | 51.6% (5.27 vs oracle 8.88 dB) |
| patch_back2back / predict | ERROR adapt_predict_init:TooFewElements | ERROR adapt_music_doa:TooFewElements |

- **Defect A: the OPT-IN adaptive loading (P9) blocks `lcmv` on small arrays.**
  It requires n_el > 2*n_comp and THROWS otherwise. P9 measured adaptive ==
  fixed at +0.00 pp on all 90 cells, so this is pure downside: no measured
  benefit anywhere, and it makes small apertures unusable.
- **Defect B: `predict` throws when n_el <= 2*n_comp** instead of degrading to
  the reactive path. MUSIC's model order is hardcoded `n_sig = 2*n_comp`.

With fixed loading the Dipole result is the ideal graceful degradation: it
matches the oracle exactly, because with one element the quiescent beam IS the
best possible. The oracle-tracking metric reports 100% -- "reached this array's
potential", and that potential is nothing. This validates the metric for the
"as close as possible to each array's potential" goal.

## O-F4. The on/off anticipatory branch essentially never fires
patchs_with_monopoles, jammer at 20 deg separation, J/N 25 dB, duty 0.5:

| toggle period | presence% | period trusted% | PRE-NULL FIRES% | lcmv | predict |
|---|---|---|---|---|---|
| 2 s | 100% | 0% | 0.0% | 37.6 | 37.6 |
| 5 s | 100% | 0% | 0.0% | 36.4 | 36.4 |
| 10 s | 81% | 58% | **1.3%** | 63.3 | 64.0 |
| 15 s | 71% | 50% | 0.1% | 75.5 | 75.6 |
| 20 s | 66% | 0% | 0.0% | 81.6 | 81.6 |
| 30 s | 60% | 0% | 0.0% | 87.8 | 87.8 |
| 40 s | 58% | 0% | 0.0% | 90.8 | 90.8 |

The entire reason P8 exists fires on <= 1.3% of steps. Two independent causes:

- **Cause A (period >= 20 s): the buffer is too short.** `buffer_len` = 1024
  steps = 51.2 s and `min_periods` = 3 requires three cycles inside the window,
  so a period above ~17 s can NEVER be trusted. Arithmetic, not tuning.
- **Cause B (period <= 5 s): presence is saturated at 100%.** The covariance
  horizon is 1/(1-lambda) = 10 steps = 0.5 s, but the OFF phase is only
  1-2.5 s, so R_hat still carries jammer energy and the eigengap never
  collapses. **The presence signal the periodogram consumes is itself low-pass
  filtered by lambda.** Constant presence -> no spectral line -> never trusted.
  (Note presence reads 81% at duty 0.5, i.e. it over-reports ON by ~31 pp.)
- Even in the 10-15 s band where detection works, `lead_steps` = 6 steps
  (0.3 s) makes the pre-null window a negligible fraction of a cycle.

**Candidate fixes (to be implemented and re-measured in O3):**
1. Derive `buffer_len` from the period range to be detectable, or analyse the
   full presence record up to a cap, instead of a fixed 1024.
2. Scale `lead_steps` with the DETECTED period rather than a fixed constant.
3. Dual-timescale presence: a short-lambda covariance for presence detection
   alongside the long-lambda one used for beamforming. This is the real fix for
   fast toggling and the deepest change of the three.
4. Adaptive-loading precondition -> degrade to fixed loading, never throw.
5. MUSIC DoF check -> fall back to the reactive path, never throw.

## O-F5. The repair works, and the mechanism numbers prove it
Opt-in `adapt.predict.onoff` block: presence from a short-memory covariance,
analysis window sized from `max_period_s`, lead scaled to the detected period
and capped at a small number of covariance horizons.

**Presence error against the true 50% duty cycle** (campaign, 137/198 cases):

| toggle period | shipped | repaired |
|---|---|---|
| 4 s | 0.455 | **0.122** |
| 10 s | 0.334 | **0.049** |
| 25 s | 0.135 | **0.020** |

Pre-null firing rate 0.000-0.010 -> 0.014-0.047. Period trust 0% -> 50-90%.

**Score, mean over cells:** predict 69.2 -> **76.7** (+7.5 pp), cells at or
above 90: 36 -> 41. `lcmv` unchanged at 66.7 (the repair only touches predict).

| period | base predict | repaired | paired mean delta | cells worse |
|---|---|---|---|---|
| 4 s | 56.9 | 66.7 | +9.8 | **15 of 46** |
| 10 s | 70.6 | 79.5 | +8.9 | 0 of 46 |
| 25 s | 80.4 | 84.1 | +3.7 | 0 of 45 |

Per array (predict, base -> repaired): patchs 87.9 -> 87.8 (already at its
ceiling), Monopoles 80.2 -> 86.3, spacing0.6 60.7 -> 70.3,
spacing0.6_disturbed3 60.4 -> 71.2 (cells passing 1 -> 4).

## O-F6. The T=4 s regression is a RELEASE-POLICY finding, not a tuning miss
All 15 regressions are at the 4 s period, all on cells the reactive path
already handled (base 78.8-89.2). Two hypotheses were tested and BOTH refuted:

- *Mis-timed pre-null?* No. Shortening `fast_lambda` 0.7 -> 0.5 -> 0.3 improves
  presence error 0.30 -> 0.15 -> 0.09 while the score stays at **exactly 80.1**.
- *Pre-null firing too often?* No. Cutting `lead_frac` 0.15 -> 0.03 halves the
  firing rate 4.7% -> 2.5% while the score stays at **exactly 80.1**.

The real cause: with presence saturated at ~100% (the pre-repair state)
`predict` NEVER reached its "jammer absent" branch, so it silently behaved as
`lcmv`. Repairing presence exposed the release policy for the first time --
on OFF it drops the null and returns to the quiescent beam. At a 4 s period the
OFF window is 2 s and turn-ons are frequent, so re-acquiring the null costs
more than the extra gain buys.

**A release gate was implemented and measured, then defaulted OFF.** Requiring
the OFF window to exceed N covariance horizons before releasing does fix the
tail (89.4 -> 90.0 instead of 89.4 -> 80.1) but destroys the wins: a hard 10 s
cell collapses 62.2 -> 33.0, because holding the null forfeits exactly the
quiescent gain that made the repair pay. Net, the un-gated form is +7.5 pp
overall against a small tail of -1 to -8.6 pp. Shipped un-gated, knob retained
and documented.

## O-F7. Process note: a live-edit hazard, caught and neutralised
`adapt_predict_update.m` was edited while both campaign arms were running, and
MATLAB reloads changed functions on the next call -- so the `onoff` arm could
have changed behaviour mid-run at ~case 137. Neutralised by defaulting the new
gate to inert and VERIFYING the new code reproduces the pre-edit numbers
exactly (80.1 / 81.0 / 82.6 / 62.2 / 80.4 on the five probe cells). The arm is
therefore homogeneous. Worth a standing rule: do not edit algorithm modules
while a campaign is in flight.
