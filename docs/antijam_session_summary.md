# Anti-jam Mode C — session conclusions (2026-09-07 … 09-12)

Consolidated outcome of three campaigns (P12b, Phase O, Phase O2) totalling
~8,700 closed-loop runs across 7 arrays. Read this before working on the
anti-jam milestone: it records what is settled, what is still open, and which
earlier statements were retracted, so none of it has to be rediscovered.

Full detail: `docs/antijam_p12b/`, `docs/antijam_phaseO/`.
Phase `Status` fields in `antijam_milestone_plan.md` remain the source of truth
for progress.

---

## 1. What the algorithms actually are

Two Mode C algorithms exist, and it is worth naming them precisely because the
customer-facing framing is easy to get wrong.

| name in code | what it does | how it finds the jammer |
|---|---|---|
| `lcmv` | Minimises output power subject to a distortionless constraint on the target (MVDR/LCMV). | **Never locates it.** The jammer is nulled as a by-product of minimising power. |
| `predict` | Same, plus MUSIC to estimate the jammer direction explicitly, then a hard null at a chosen angle. | Explicit direction estimate. |

**Important for any presentation:** `lcmv` is *not* a search over SINR. It is a
closed-form solution that is mathematically equivalent to maximising SINR, but
nothing iterates or "tries" weights. The only genuinely direct-search method in
the tree is `spsa` (Mode S), which probes a scalar SINR feedback channel and is
slow (~1000 probes on a real array). Describing `lcmv` as "directly optimising
SINR" is defensible only with that clarification.

---

## 2. Settled conclusions

**2.1 The drifting-jammer failure was angle lag, of exactly one covariance
horizon.** A null steered from ground truth but delayed 10 steps reproduces
`lcmv`'s score; 10 = 1/(1−λ) at λ = 0.90. Steering at the true current angle
gives 94.8% against `lcmv`'s 44.8%. A constant-velocity DoA predictor
(`adapt_cv_*`) closes it: DRIFT 55.2 → 90.3 mean, cells passing 2 → 11 of 15,
with STATIC and on/off runs **bit-identical** (max |ΔSINR| < 1e-12) because the
branch is gated on estimated angular speed. Re-tuning λ cannot substitute — no
value makes the failing cells pass.

**2.2 The on/off predictor was never running.** Its anticipatory branch fired on
≤ 1.3% of steps, and not at all outside a 10–15 s band. Two arithmetic causes:
the analysis window could not hold `min_periods` cycles above ~17 s, and below
~5 s presence saturated at 100% because the detector read the same
λ-smoothed covariance the beamformer uses. Fix: a second short-memory covariance
used **only** for detection. Presence error against the true duty fell
0.46/0.32/0.14 → 0.13/0.06/0.03 at 4/10/25 s.

**2.3 The release-policy question is CLOSED.** Taking the better of {no repair,
binary release} per cell — an oracle no realisable rule can beat — is worth
**+0.9 mean and zero additional passing cells**. A graded release was built and
measured worse (76.1 vs 80.6, worse on 45 of 90). Do not spend further effort on
*when to drop the null*.

**2.4 Calibration is the binding constraint on fielding, not any algorithm
choice.** With per-element phase error the usable region is below **~1° RMS**;
by 3° every method is at the floor, and the p10–p90 spread at 2° is 1.1–76.8, so
a calibration statement needs a distribution rather than a point. This was
structurally unmeasurable until an opt-in assumed-patterns path was added to
`closed_loop_run`.

**2.5 `guard_deg = 5` is wrong on every array.** Derived per (array, target)
from the beamwidth it spans **16°–72.5°**. Any earlier result whose jammer sat
between the true and assumed guard was measuring a main-beam jammer, which the
milestone declares out of scope. `kpi_array_profile` derives it.

**2.6 The P9 data-driven loading buys nothing.** Adaptive minus fixed diagonal
loading measured **+0.00 pp on all 90 cells**, in both comparisons — identically
zero every cell. It also used to *throw* on small apertures, taking the reactive
tracker down with it. Pure downside.

**2.7 Generality: every array now runs, and infeasible ones say so.** A
1-element and a 2-element array previously crashed the stack. They now degrade
with a warning. A geometry preflight refuses 11 of 35 (array, target) pairs with
a stated reason. The 1-element array scores 100.0 (5.92 dB against an oracle of
5.92 dB) — with no spatial degrees of freedom the quiescent beam *is* the
optimum, and the metric correctly reports the array reached its potential.

**2.8 Planar arrays have an exact up/down ambiguity.** `e(θ,φ) ≡ e(180−θ,φ)` to
2.8e-6 on `ManyDipoles`. **Benign for nulling** (identical steering vectors, so a
null at the mirror is the null at truth) but **fatal for tracking**. Measured
modulo the fold, that array's DoA error is 1.5° median, not 139°. Any DoA-error
KPI there is meaningless unless computed modulo the fold.

**2.9 Availability is not a usable headline metric.** A 35-point improvement in
oracle-tracking score moved availability by **+0.05 pp** (98.18 → 98.23). It
saturates on easy cells and cannot rank algorithms. Use the oracle-tracking
score; keep availability as a diagnostic.

**2.10 Score and potential must always be reported together.** A 2-element array
scores 99 because its achievable SINR is 11 dB; a 16-element array scores 72
against an achievable 31 dB. A high closeness score is *not* good absolute
performance.

---

## 3. Retractions made this session

- **Phase O's "the on/off repair regresses ManyDipoles, do not enable"** was
  **wrong** — a run-length artifact. Period learning needs 3 cycles and the runs
  were 4 cycles long. At 8 cycles the same array reads 73.6 → 77.8 and the repair
  improves or ties on every array. Now enabled by default. The risk had already
  been written down as a known limit before the recommendation was issued.
- **The 2026-09-06 front/back-ambiguity retraction** measured the wrong angle
  pair (antipodal (180−θ, φ+180) instead of the mirror (180−θ, φ)). The
  ambiguity is real.

---

## 4. Still open

- Calibration robustness (robust MVDR / uncertainty-tied loading). Highest-value
  remaining algorithm work.
- MUSIC model order hardcoded at `2·n_comp`; a rank-aware estimate would remove a
  class of array-specific failures.
- MUSIC cost is 42–126× the reactive path — blocking for any real-time target.
- `guard_deg` derived in the harness but still a global constant in `config.yaml`.
- Not qualified above ~10°/s drift.
- Out of scope but plausible in the field: multiple simultaneous jammers,
  polarization mismatch.

---

## 5. Method lessons that generalise

- **When a score refuses to move as you sweep the supposed cause, the cause is
  elsewhere.** Two confident explanations of a regression were falsified this way
  — the score sat at exactly 80.1 under two different knob sweeps.
- **Weight a tuning set by the population's own incidence.** A policy tuned on
  five hand-picked conflict cells (60% problem cases) lost at scale, where
  problems are 11 of 90 (12%).
- **Run length can decide the conclusion.** See §3.
- **Do not edit algorithm modules while a campaign is in flight** — MATLAB
  reloads changed functions on the next call and can split an arm in two.
- **`NaN <= 0` is false.** A guard written that way let a run complete while
  silently returning garbage (a cell read 12.8 instead of 88.2). Test finiteness
  explicitly.
