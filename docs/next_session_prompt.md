# Prompt for the "clear rewrite" session

Copy the block below as the opening message of a new session.

---

I need a clean, readable reimplementation of the anti-jam beamforming work in a
**new folder**, `MATLAB/antijam_clear/`, kept completely separate from the
existing `MATLAB/antijam_utils/` tree. Do not modify anything outside the new
folder except where I say so below.

**Why.** I have a customer meeting in a few days. I do not need to hand over
code at that meeting — I need to *present what has been done and explain how it
works*. The existing implementation performs well but I cannot walk someone
through it, and I will not present something I cannot explain. **I would rather
have a simpler algorithm I fully understand than a better one I don't.** Poorer
numerical results are an acceptable trade; unexplainable code is not.

Read `docs/antijam_session_summary.md` first — it records what three campaigns
established, including which questions are already closed, so you do not
rediscover any of it.

## What to build

Two algorithms, deliberately separate, each in its own file, each doing one
thing in a way I can follow:

1. **`DetectAndNull`** — the explicit approach. Estimate *where* the jammer is
   and *what it is doing* (steady / on-off / drifting linearly), then shape the
   radiation pattern accordingly. The logic should be legible as: detect, then
   decide, then place the null.

2. **`MaxSinrWeights`** — the implicit approach. Solve directly for the weights
   that maximise output SINR, without ever locating the jammer. This is the
   classical MVDR/LCMV solution and it is closed-form, not a search.

The top-level script should read like the structure below, with each step a
named call whose purpose is obvious from its name:

```matlab
%% Load the array and the scenario
[element_patterns, theta_deg, phi_deg] = LoadArray(array_folder);
signal_angle = ...;

%% Approach 1 - detect the jammer, then null it
jammer_state  = DetectJammer(snapshots, element_patterns, theta_deg, phi_deg);
weights_1     = DetectAndNull(element_patterns, signal_angle, jammer_state);

%% Approach 2 - optimise SINR directly, no detection
weights_2     = MaxSinrWeights(snapshots, element_patterns, signal_angle);

%% Compare
CompareApproaches(weights_1, weights_2, ...);
```

Adjust the names if something reads better, but keep that shape: a person should
be able to read the main script top to bottom and understand the whole method
without opening any other file.

## Rules

- **Explain every logical step as you go.** Before writing each function, tell me
  in plain language what it does, why it is needed, and what would go wrong
  without it. I want to understand the reasoning, not just receive code. Pause
  and let me confirm I follow before moving on — I would rather this take longer
  and land.
- **Simplicity beats performance.** If a simpler formulation costs a few dB, take
  it and say what it cost.
- **No machinery.** No opt-in config blocks, no campaign runners, no arm/variant
  sweeps, no tuning knobs beyond what is genuinely needed. If a parameter exists,
  I should be able to say what it does and why it has that value.
- **Reuse, do not re-derive:** the CST parser and the data in `data/` are fine as
  they are. Use `matlab_utils/` for loading. Do not depend on `antijam_utils/`.
- MATLAB, must run on **R2020a** (no `clim`, no `arguments` blocks).
- Follow `docs/STYLE.md`, and CLAUDE.md's rules — especially: ask before assuming,
  and describe an approach before coding it.

## Things that are already known — respect these, don't rediscover

These cost real campaign time to establish. Details in the summary document.

- The guard sector around the target must be **derived from the array's own
  beamwidth** per (array, target). The global 5° in `config.yaml` is wrong on
  every array (true range 16–72.5°).
- Diagonal loading is **required**, not optional: the snapshots contain the
  desired signal, so without it a steering error makes the beamformer null your
  own signal.
- Some arrays (e.g. `ManyDipoles`) are **exactly ambiguous between θ and 180−θ**.
  Harmless when nulling, ruinous when tracking a direction over time.
- `Dipole` (1 element) and `patch_back2back` (2 elements) cannot null. Handle
  them gracefully — they are a good demonstration that the method knows its own
  limits.
- Score against **what that array can achieve**, not an absolute target, and
  always report the achievable figure next to the score.

## Deliverables

1. The code in `MATLAB/antijam_clear/`, with a runnable demo script.
2. A walkthrough document explaining each step, at the level of the existing
   `docs/antijam_phaseO/explain_technical.html` but covering only this simpler
   implementation.
3. Presentation material I can use in the meeting: the logic of both approaches,
   a worked example on one array, and an honest statement of what this simpler
   version gives up against the existing one.

## Done means

I can read the main script and explain, out loud and without notes, what every
step does and why it is there.

Start by proposing the file layout and the list of functions with one line each
on what they do. Do not write code until I have agreed to that list.
