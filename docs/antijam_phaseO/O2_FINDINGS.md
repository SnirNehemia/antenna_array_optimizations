# Phase O2 — graded release. Running findings.

## O2-F1. The binary release could not work, and the graded one needs a learned period
Phase O left a release policy that was all-or-nothing: the instant presence went
false, the null was dropped and the beam returned to quiescent. Measured then:
it wins where the gap is long and loses where it is short, and NO fixed
threshold can serve both, because any threshold holds through the first N steps
of EVERY OFF window -- including the long ones where releasing is right.

**The replacement relaxes the null CONTINUOUSLY via diagonal loading.** Loading
is the physical knob that trades null depth for main-beam gain, so raising it
toward the scale of R_hat's own energy progressively swamps the jammer
eigenvalue and returns the beam to quiescent. Chosen over blending two weight
vectors, which are defined only up to a phase and can cancel.

## O2-F2. A fixed ramp still trades; an OFF-WINDOW-ADAPTIVE ramp does not
Probe over 5 conflict cells (3 fast-toggle regressions, 2 long-period wins):

**4 cycles** (period rarely learned):
| config (start/ramp/off_frac) | mean | regress cells | win cells |
|---|---|---|---|
| 0/0/0 (un-graded) | 69.5 | 65.9 | 73.9 |
| 3/5/0 fixed | 74.1 | 83.4 | 62.5 |
| 2/4/0 fixed | 72.8 | 79.2 | 64.9 |

A fixed ramp buys the regressions back and pays for it in the wins -- better
balanced than binary, but still a trade.

**12 cycles** (period learned; the adaptive ramp can size itself):
| config | ManyD (90,30) T4 | ManyD (45,150) T4 | patchs T4 | sp0.6 T10 | dist3 T10 | mean |
|---|---|---|---|---|---|---|
| base (no repair) | 97.8 | 75.5 | 87.7 | 28.4 | 60.2 | 69.9 |
| un-graded | 79.0 | 61.3 | 86.6 | **55.9** | **82.6** | 73.1 |
| fixed 3/5 | 97.8 | 90.1 | 91.7 | 36.5 | 64.5 | 76.1 |
| **adaptive 2/2/0.30** | 96.5 | 86.9 | 90.7 | 42.5 | 69.6 | **77.2** |

The adaptive ramp beats the un-graded repair AND the fixed ramp, and no longer
falls below base on the fast-toggle cells. It only works once the period is
learned, which needs more than 4 cycles -- confirming O-F4 (that the previous
campaign's 4-cycle runs were an artifact) as a real limitation rather than a
caveat.

## O2-F3. A NaN guard that a gate now enforces
With the onoff block ABSENT the release fields are NaN, and `NaN <= 0` is false
-- so the ramp maths ran on NaN, the diagonal loading went NaN, and the
beamformer solve turned singular. The run still COMPLETED, silently returning
garbage: a 25 s cell read 12.8 instead of 88.2. Caught only because the probe
printed the un-repaired baseline alongside.

Two lessons taken: the guard is now explicit about finiteness, and
`test_antijam_onoff` asserts that an un-repaired run produces finite SINR and
finite weights at every step -- which is what would have caught it.

## O2-F4. Campaign re-run
90 testable cases x 3 seeds x 3 algorithms x 3 arms {base, onoff, graded}.
**8 cycles per run** rather than 4, so the period is learned and the adaptive
ramp can engage. Grid trimmed to 3 targets x 2 separations to pay for the
doubled run length; MUSIC stride raised to 4 on the 1 deg arrays (DoA RMSE
1.00 -> 1.45 deg, acceptable where the question is presence rather than angle).

## O2-F5. VERDICT: the graded release is WORSE at scale. Negative result.
90 cases x 3 seeds x 3 algorithms x 3 arms, 8 cycles, 2,430 runs, 0 failed
case-seeds:

| arm (predict) | mean score | cells >= 90 of 90 |
|---|---|---|
| base (no repair) | 73.5 | 24 |
| **graded release** | 76.1 | 29 |
| **binary release (Phase O repair)** | **80.6** | **32** |

Paired per cell, on predict:
- binary vs base: **+7.10** mean, 69 better / 11 worse, worst -25.5
- graded vs base: +2.63 mean, 66 better / 10 worse, worst -14.0
- **graded vs binary: -4.46 mean, 9 better / 45 WORSE**, worst -31.5

Binary wins at EVERY toggle period (4 s: 69.0 vs 62.7; 10 s: 83.1 vs 78.1;
25 s: 89.8 vs 87.6) and on every array except a 0.9 pp edge on ManyDipoles.

**Why the probe misled me.** The 5-cell tuning set was hand-picked as "conflict
cells" -- 3 fast-toggle regressions against 2 long-period wins. That is a 60/40
split of a population that is really 11/79. The graded ramp's cost falls on the
majority (partially holding the null through OFF, forfeiting gain the oracle
already has) and its benefit on the minority. Tuning on a curated conflict set
over-weighted exactly the cells that motivated the change. A cheap guard would
have been to weight the probe set by the campaign's own incidence.

## O2-F6. There is almost nothing left to win from release policy -- question closed
Taking the better of {base, binary} PER CELL -- an oracle no real rule could
beat -- scores **81.5 mean / 32 passing** against binary's **80.6 / 32**.

So any switching, hybrid or confidence-weighted release policy is bounded at
**+0.9 mean and zero additional passing cells**. The graded release was worth
building to establish that number, and the number says stop: further effort
belongs elsewhere (detection quality, MUSIC model order, calibration), not in
deciding when to drop the null.

## O2-F7. RETRACTION: Phase O's "regresses ManyDipoles" was a run-length artifact
Phase O reported the binary repair regressing ManyDipoles 83.0 -> 79.7 and
recommended against enabling it. At 8 cycles instead of 4 that reverses:
**73.6 -> 77.8 (+4.2)**, and the binary repair now improves or ties on EVERY
array:

| array | potential | lcmv | base | graded | binary |
|---|---|---|---|---|---|
| spacing0.6 | 30.7 dB | 58.5 | 61.6 | 65.6 | **72.1** |
| spacing0.6_disturbed3 | 29.3 dB | 56.7 | 60.7 | 65.1 | **75.5** |
| Monopoles | 19.7 dB | 77.5 | 79.5 | 78.2 | **84.5** |
| patchs_with_monopoles | 25.9 dB | 90.4 | 90.8 | 92.0 | **92.0** |
| ManyDipoles | 22.2 dB | 73.6 | 73.6 | **78.7** | 77.8 |
| patch_back2back | 11.0 dB | 98.9 | 98.9 | 98.9 | 98.9 |

O-F4 predicted this: period learning needs 3 cycles, so a 4-cycle run spends
most of itself unlearned and understates the anticipatory path. The prediction
was right and the earlier recommendation was wrong.

**Recommendation now REVERSED: enable the binary on/off repair by default.**
Recovery after turn-on also improves (14.7 -> 9.7 steps at 10 s, 25.1 -> 16.5
at 25 s). Keep the graded knob available for a deployment known to face fast
toggling on a coarse-grid array, where it is the better choice (ManyDipoles at
a 4 s period: binary -25.5, graded +0.2).
