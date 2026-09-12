# Anti-jam milestone — rewrite brief

**Purpose.** The anti-jam milestone (P0–P12b, Phase O, O2) is functionally complete
and its conclusions are settled. It grew phase by phase, and it shows: 7,665 lines in
`MATLAB/antijam_utils/`, of which the algorithm proper is ~2,200 and the figure code is
~2,950. This document is the input to a **rewrite** that keeps the physics and the
measured conclusions but rebuilds the structure.

It is written for someone starting cold. It says what the milestone is actually about,
what earned its place and what did not, which problems are *closed* (so the rewrite
does not relitigate them), and which mistakes cost the most time — because most of them
are structural and a rewrite can design them out.

**The rewrite's goal is explainability, not capability.** `docs/next_session_prompt.md`
is the brief for that session: a clean reimplementation in `MATLAB/antijam_clear/`, two
deliberately separate approaches (detect-then-null, and the closed-form max-SINR
solution), readable top to bottom, for a customer meeting where the work has to be
*explained* rather than handed over. **A simpler algorithm that is fully understood beats
a better one that is not**, and a few dB is an acceptable price. Read this document for
the physics and the settled facts; read the prompt for the shape of the deliverable.

Read alongside: `antijam_milestone_plan.md` (phase-by-phase record of what was built and
why), `docs/notes.md` (session log; the `[Pn]`/`[O]` prefixes index it),
`docs/antijam_phaseO/onoff_report.html` (the on/off results in narrative form),
`docs/antijam_p12b/modec_deck.html` (the whole milestone in 22 slides).

---

## 1. The irreducible problem

One jammer. Unknown direction. 2-D angles (θ, φ). Real CST element patterns. At each
step, choose a complex weight per element.

```
SINR(w) = σs²|wᴴe_s|² / (σj²|wᴴe_j|² + σn²‖w‖²)
```

Two observation regimes, and an algorithm belongs to exactly one:

- **Mode C** — per-element complex snapshots. The richer regime, and the one a real
  digital-beamforming receiver offers. Everything interesting is here.
- **Mode S** — a scalar SINR reading only. SPSA and the bandit live here.

**The one inviolable rule**: an algorithm may consume only the `obs` struct for its
declared mode. The true jammer angle goes to the oracle and nowhere else. The oracle is
not deliverable; it exists to bound the others.

---

## 2. What earned its place

Keep these. Each is load-bearing and each is backed by a measurement, not a preference.

### 2.1 The closed-form beamformer

```
w = R⁻¹e(θs) / (e(θs)ᴴ R⁻¹ e(θs))     with diagonal loading
```

Everything else in the stack is either **feeding R** or **aiming e**. If the rewrite has
a spine, this is it. Note it is the **MPDR** form — the wanted signal is inside R. That
is convenient in practice and is also the cause of the largest open problem (§5.1).

### 2.2 One time constant, `forgetting_lambda = 0.90`

`R̂(t) = λR̂(t−1) + (1−λ)xxᴴ`, applied **once per step to the batch-averaged snapshot
covariance** — not K sequential per-column recursions, which is a bug the original code
had and which made λ = 0.98 look correct.

The horizon `1/(1−λ) = 10` steps explains an implausible fraction of everything observed
across three campaigns:

- the drifting-jammer failure (the null aims 10 steps into the past);
- the derived lead of 14 steps = 10 + 4 (§2.4);
- why the on/off presence detector was blind (it read a covariance that smooths the very
  edge it was trying to find);
- why a null usefully survives a jammer switching off.

**Design rule: if λ moves, the derived constants move with it.** A rewrite should make
that dependency explicit in code rather than in a comment.

### 2.3 The oracle-tracking score

```
track_score = 100 · mean( (SINR_oracle − SINR) ≤ 3 dB )
```

The fraction of a run spent within 3 dB of the perfect-knowledge beamformer, **normalised
per cell**, so a 6-element patch array and a 20-element dipole array are each scored
against their own potential. 90 is the pass mark (a judgement, stated so it can be
re-read against a different one).

This replaced **availability**, which saturates on easy cells: it moved +0.05 pp for a
change worth 35 points on this metric. Availability cannot rank algorithms. Do not
bring it back as a headline.

A pleasant property worth preserving: an array that can do nothing scores 100 for doing
nothing. The 1-element array's quiescent beam *is* its optimum (5.92 dB against an oracle
of 5.92 dB), and the metric says so rather than marking it failed.

**Always report the score and the potential together.** A 2-element array scores 99
because its achievable SINR is 11 dB; a 16-element array scores 72 against an achievable
31 dB. A high closeness score is *not* good absolute performance, and a table of scores
alone will mislead every reader who sees it.

### 2.4 The CV-Kalman DoA predictor

A constant-velocity Kalman filter on `[θ, θ̇, φ, φ̇]`, fed by the MUSIC angle the stack
already computes, extrapolating `lead = 14` steps. Worth **+35 points** on drift
(55.0 → 90.1) and provably inert everywhere else — gated on estimated angular speed, so
STATIC and WINDOW runs are bit-identical with it on or off.

Two details that must survive the rewrite:

- **The lead is derived, not fitted.** 10 (covariance horizon) + 4 (one-step application
  delay plus the `weight_smoothing_mu` first-order lag). A tuned lead would have been a
  much weaker result.
- **θ-mirror folding.** Planar arrays have an exact up/down ambiguity:
  `e(θ,φ) ≡ e(180−θ,φ)` to 2.8e-6 on `ManyDipoles`. It is **benign for nulling** — the
  steering vectors are identical, so a null at the mirror *is* the null at truth — and
  **fatal for tracking**. The filter probes for the degeneracy
  (`mirror_coherence > 0.99`) and folds rather than chasing a phantom. Consequence for
  measurement: **any DoA-error KPI on such an array is meaningless unless computed modulo
  the fold.** Measured modulo the fold, that array's DoA error is 1.5° median, not 139°.

### 2.5 The fast presence covariance

A **second** EMA at a shorter λ (`fast_lambda = 0.5`), used **only** for presence
detection and never to form weights. One extra line of state that took the anticipatory
branch from firing on ≤1.3% of steps to 1.4–4.7%, and presence error against the true
duty cycle from 0.46 / 0.32 / 0.14 to 0.13 / 0.06 / 0.03 at toggle periods of 4 / 10 /
25 s.

The principle generalises: **a detector slower than the beamformer cannot resolve edges
the beamformer already smooths.** Any future detector gets its own memory.

### 2.6 Graceful degradation, everywhere

No module may throw because an array is small or degenerate. A 1-element and a 2-element
array used to crash the entire stack. Now:

- `adapt_music_doa` returns `feasible = false` instead of throwing when `n_sig >= n_el`;
- the loading precondition warns and degrades instead of throwing;
- a geometry preflight refuses 11 of 35 (array, target) pairs as **not real tests**,
  each with a stated reason, rather than scoring them as failures.

### 2.7 `kpi_array_profile` and the derived guard

Per (array, target): directivity, HPBW in θ and φ, **derived `guard_deg`** (half the
wider HPBW), mirror coherence, MUSIC feasibility. This is what makes cells comparable
across arrays, and it supplies the feasibility check above.

**Steering coherence** between the signal and jammer directions is the array-independent
difficulty coordinate. Use it as the x-axis whenever comparing arrays.

---

## 3. What did not earn its place

### 3.1 Closed questions — do not reimplement

- **Graded release** (relaxing the null continuously instead of dropping it). Built,
  measured at scale, rejected: 76.1 against the binary policy's 80.6, worse on 45 of 90
  cells. More importantly, **best-of-{none, binary} per cell — an oracle no real
  switching rule can beat — scores 81.5 against binary's 80.6.** The entire family of
  switching, hybrid and confidence-weighted release policies is bounded at **+0.9 mean
  and zero extra passing cells.** That bound is the useful output of having built it.
  The code ships present and disabled; the rewrite can drop it to a documented note.

- **Data-driven adaptive loading** (`loading_factor_db`). Measured at **+0.00 pp on all
  90 cells**, and it was the direct cause of two arrays being unable to run at all.
  Delete it. (Note: an earlier session wrongly blamed it for a calibration collapse and
  that claim was retracted — it is not harmful, it is merely worthless.)

- **λ re-tuning as a fix for drift.** Swept; there is no value at which the failing
  cells pass. The deficit is angle lag, not smoothing.

### 3.2 Dead weight relative to a Mode C rewrite

- **Mode S** — `adapt_spsa_*`, `agent_bandit_*`, `agent_codebook_build` (367 + ~200
  lines). A different regime sharing a folder. Decide explicitly whether it is in scope;
  if not, it should not be in the same directory.
- **The frequency / notch layer** — `adapt_freq_*`, `sim_notch_response` (~370 lines).
  Entirely inert unless `sim.fs_hz` turns the waveform layer on. Working and tested, but
  orthogonal.
- **Figure code** — 12 `plot_*` and 4 `save_*` modules, **2,950 lines, 38% of the tree**,
  grown one-per-phase with no shared style layer. This is the single largest
  consolidation opportunity and the lowest-risk one.

### 3.3 `config.yaml` is a changelog wearing a config's clothes

The `adapt:` block is a few dozen keys buried in several hundred lines of "tried X on
date Y, got Z, kept W." That history is genuinely valuable — it is why the numbers are
what they are — but it belongs in `notes.md` or the plan, with the config carrying only
what a reader needs to *set* the key.

---

## 4. Structural problems to design out

1. **`adapt_predict_update.m` is 485 lines** and performs: two covariance updates, MUSIC,
   presence detection, period estimation, lead computation, CV prediction, release
   policy, weight formation, and smoothing. It should be a short pipeline of named
   stages, each independently testable.

2. **The opt-in config-block pattern is the wrong inheritance.** "Absent block ⇒
   byte-identical legacy behaviour" was exactly right for de-risking incremental phases
   against a frozen gate suite. It is exactly wrong as a starting structure: it is why
   fields exist in a NaN state at all, and it produced the worst bug of the project
   (§5.2, trap 1). A rewrite has no legacy to be identical to. **Every field should be
   constructed valid or not constructed.**

3. **MUSIC's model order is hardcoded** at `n_sig = 2·n_comp`. An assumption, not a
   measurement. It makes MUSIC impossible on small dual-polarization apertures and
   mis-reads presence on arrays whose second component is weak but real. Estimate the
   order from the eigenvalue profile and a whole class of array-specific failure
   disappears.

4. **`guard_deg` ships as a global 5°** while the true per-(array, target) value spans
   **16°–72.5°**. `kpi_array_profile` already derives it. Any result whose jammer sat
   between the assumed and the true guard was measuring a jammer inside the main beam —
   out of scope by the milestone's own definition.

5. **Gate suites are self-contained**, each with its own hardcoded config, so
   `config.yaml` is untested surface. The isolation is worth keeping; the gap should be
   closed by a separate contract test over the shipped config, not by coupling the
   gates to it.

---

## 5. What is genuinely open

Ordered by value. This is the rewrite's real target list.

### 5.1 Steering mismatch at strong signal — the largest named number

Sweeping signal and jammer amplitude independently:

| array | of the 20 dB the oracle gains, we capture | lost |
|---|---|---|
| spacing0.6 | 11.9 dB | **8.1** |
| patchs_with_monopoles | 17.3 dB | 2.7 |
| ManyDipoles | 19.7 dB | 0.3 |
| Monopoles | 19.8 dB | 0.2 |

Raising the **jammer** 20 dB costs at most 2.1 dB — the nulling works. Raising the
**wanted signal** 20 dB costs up to 8.1 dB. The mechanism is MPDR self-cancellation: the
wanted signal is inside R, so under steering mismatch the solver spends degrees of
freedom cancelling it, and the stronger it is the more there is to cancel. The oracle is
immune because it is handed the true steering vector.

Two things follow, and both matter:

- **Every campaign number in the milestone is a weak-signal number** (σs = 0 dB) and so
  flatters the stack relative to a strong-signal deployment.
- **8.1 dB against the +0.9 dB a perfect release policy was worth.** This is where the
  effort belongs.

Standard remedies, cheap to try in this order: diagonal loading swept against *mismatch*
rather than against null depth; a derivative (point-and-null) constraint widening the
protected region around θs; robust Capon with an explicit uncertainty ellipsoid if those
fall short. **Measure it with the amplitude grid, not the campaign** — the campaign runs
at σs = 0 dB and structurally cannot see the effect.

**Not established**: why the loss ranges 8.1 → 0.2 dB across four arrays. The worst case
being the highest-directivity array fits the mismatch story (a narrow beam pays more for
the same angular error), but four arrays at one geometry each is an observation, not a
controlled sweep. A mismatch × σs sweep would settle it.

### 5.2 Calibration tolerance

Independently measured: the stack is usable below **~1° RMS** per-element phase error and
at the floor by 3°. The spread is enormous — at 2° the p10–p90 range is 1.1–76.8, so one
array build could be fine and the next unusable. Unchanged by the predictor, which keeps
its 2–4× advantage throughout. This is an exposure of the whole approach, and §5.1 is the
same weakness measured from the other direction.

### 5.3 Smaller, well-specified

- Derive `guard_deg` per (array, target); retire the constant. Cheap.
- Rank-aware MUSIC model order. Moderate.
- **MUSIC costs 42–126× the reactive path.** Blocking for any real-time target, and
  unaddressed. If the rewrite has a compute budget, this is the number it collides with —
  and it is an argument for making the DoA stage optional by construction rather than
  always-on.
- Not qualified above **~10°/s** drift — the predictor extends the envelope, it does not
  remove the limit.
- Multiple simultaneous jammers; polarization mismatch; real-time compute budget.

---

## 6. Traps

These cost real time. Most are structural, and a rewrite is the moment to design them out.

1. **NaN does not compare the way you want, and a bad run can look like a good one.**
   With the `onoff` block absent its release fields were NaN; `NaN <= 0` is false, a
   guard fell through, loading went NaN, the linear solve went singular — and the run
   **completed silently with garbage** (a cell read 12.8 dB instead of 88.2). Nothing
   errored. **Every gate must assert finite SINR and finite weights, not merely that the
   run completed.** A test that only checks "it ran" would have passed this.

2. **Run length is a confound, and it shipped a wrong recommendation.** Period detection
   needs 3 cycles; the Phase O campaign ran 4, so ~75% of each run was spent before a
   period existed. That produced a regression that became the recommendation *not* to
   enable the repair. Re-run at 8 cycles it improves or ties on every array, and the
   recommendation was reversed. **State the convergence requirement of every estimator
   and make the run length a derived quantity, not a constant.**

3. **Never edit algorithm modules while a campaign is in flight.** MATLAB reloads changed
   functions on the next call, so an edit at case 137 can split an arm in two. This
   happened; it was neutralised by defaulting the new gate inert and verifying the new
   code reproduced pre-edit numbers exactly — but the rule is simply don't.

4. **Metric invariance is a diagnostic, not a nuisance.** Twice, a hypothesis was
   falsified by finding the score *invariant* to the knob being swept — which is what
   pointed at the real cause. When a metric refuses to move, the cause is elsewhere.

5. **Weight a tuning set by the population's incidence before trusting it.** The graded
   release was tuned on five hand-picked conflict cells — a 60/40 split of a population
   that is really 11/79. The benefit fell on the minority and the cost on the majority,
   and the tuning looked excellent throughout.

6. **Two MVDR solutions are defined only up to a phase and can cancel.** This is why the
   release policy blends *loading* rather than weights. Any future interpolation between
   two weight vectors must phase-align first.

7. **MATLAB R2020a.** Use `caxis`, not `clim` (R2022a). No `defaultAxesTitleColor`, no
   `WordWrap`, no `jsonencode` PrettyPrint. The verify harness is `-batch`.

8. **The figure theme leaks.** The session theme overrides figure `Color` for `print`
   regardless of the property; force colours explicitly (axes, legend text, colorbar
   label and ticks) or use `exportgraphics` with an explicit background.

9. **Sequential `matlab -batch` renders exhaust memory** even when each fits alone —
   Windows has not reclaimed the previous process's pages when the next starts. Pause
   between them. Separately, decimate pattern grids for *display* (a 181×360 `pcolor`
   × 8 panels is ~500k patches per frame and `getframe` dies); compute on the full grid.

10. **Verify rendered output by reading it back from the artifact**, not from the render
    buffer — a one-frame GIF check produced a monochrome artifact that looked like a real
    bug, and a heredoc patch once reported success while writing nothing. For the PPTX,
    exporting every slide to PNG caught collisions invisible in the XML.

---

## 6a. Things that were believed and turned out to be wrong

Recorded so the rewrite does not re-derive the wrong answer and then re-retract it.

- **"The on/off repair regresses `ManyDipoles`, do not enable it."** Wrong — a run-length
  artifact (trap 2). At 8 cycles that array reads 73.6 → 77.8 and the repair improves or
  ties everywhere. The repair is now on by default. Note the risk had been written down
  as a known limit *before* the recommendation was issued, and was still not applied to
  it.
- **"The front/back ambiguity is not real."** That retraction measured the wrong angle
  pair — the antipodal (180−θ, φ+180) instead of the mirror (180−θ, φ). The ambiguity is
  real (§2.4).
- **"P9 adaptive loading caused the calibration collapse."** Wrong; it was a different
  random error draw. Adaptive loading is worthless (§3.1), not harmful — a distinction
  worth keeping straight if anyone revisits it.

---

## 7. What this means for the clear rewrite

`docs/next_session_prompt.md` sets the shape; this section says which of the facts above
bind it and which do not.

**Must survive, because they are physics or they are cheap and load-bearing:**

- diagonal loading — **not optional** (§2.1): the snapshots contain the desired signal,
  so without loading a steering error makes the beamformer null the wanted signal;
- the guard sector derived from the array's own beamwidth, never the global 5° (§2.7);
- graceful handling of the 1- and 2-element arrays (§2.6) — they are a good demonstration
  that the method knows its own limits;
- scoring against what that array can achieve, with the achievable figure always printed
  beside the score (§2.3);
- awareness of the θ/180−θ ambiguity wherever a direction is tracked over time (§2.4).

**May be dropped, and dropping them should be stated rather than hidden:**

- the CV-Kalman predictor (§2.4) — it is the single largest complexity in the stack and
  it buys one scenario. Leaving it out costs ~35 points on drift and nothing else. That
  is a legitimate trade for a readable implementation, provided the cost is named;
- the fast presence covariance and the whole on/off anticipation path (§2.5) — same
  argument, and simpler still: without it the reactive solution is what you get;
- every campaign runner, arm sweep and opt-in config block (§4.2);
- the entire Mode S side and the frequency/notch layer (§3.2).

**Must not be reintroduced under any framing:** graded or confidence-weighted release
policies (closed, bounded at +0.9) and data-driven adaptive loading (+0.00 pp, and it
crashed two arrays). See §3.1.

**On the gate suite** (27 files, 3,535 lines, 59 gates): do not port its structure. Port
the *assertions that encode physics* — null depth, covariance convergence, the oracle
bound — and add finite-value assertions everywhere, per trap 1. A readable
implementation with three honest tests is worth more here than fifty inherited ones.

---

## 8. Reference numbers

**These are reference points, not acceptance criteria.** The clear rewrite is expected to
be worse, and that is the agreed trade — but it should be worse *by an amount someone can
state*. Measure against these and report the difference honestly; a simpler method that
gives up 5 dB and says so is a better outcome than one that quietly gives up 15.

If a rewrite *were* aiming at parity, these are what it would have to come back with:

| measurement | value |
|---|---|
| drift cells at the pass mark, before → after CV | 2 → 11 of 15 |
| drift score, lcmv / predict / predict+CV | 55.2 / 55.0 / 90.1 |
| on/off mean score, no repair → binary repair | 73.5 → 80.6 |
| on/off cells at the pass mark | 24 → 32 of 90 |
| graded release (rejected) | 76.1 |
| best-of-{none, binary} oracle bound | 81.5 |
| per-array on/off score (8 cycles) | spacing0.6 80.1 · ManyDipoles 81.0 · Monopoles 82.6 · patchs 80.4 |
| presence error, before → after | 0.46 → 0.13 |
| calibration usable limit | ~1° RMS per-element phase |
| strong-signal capture of 20 dB | 11.9 / 17.3 / 19.7 / 19.8 dB (§5.1) |

Seed-to-seed spread of the score was median **0.65 pp** (p90 4.27, max 10.46) over 198
cases, so three seeds is adequate for differences of this size. One unrelated repository
test (`test_metrics`) fails on a pre-existing MATLAB-versus-Python discrepancy and was
failing before any of this work began.

---

*Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone. Written
2026-09-12, at the close of Phase O2.*
