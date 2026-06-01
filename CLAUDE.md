# CLAUDE.md

Antenna array pattern optimization tool. Gradient descent on complex element
weights to shape array radiation patterns (peaks + nulls) from CST Studio
far-field exports. Python now; will port to MATLAB once results are validated.

## Authoritative References

Read these before working on the relevant area. Do not duplicate their content here.

- **`docs/pipeline.md`** — data flow, file formats, array model, cost function math
- **`docs/STYLE.md`** — coding style, naming, docstring format, MATLAB porting flags
- **`docs/notes.md`** — physics assumptions, open questions, session log
- **`config.yaml`** — user-facing configuration schema

## Project Map

```
src/io/         CST .txt parser
src/cost/       Objective function (peak/null directives → scalar J)
src/optimize/   L-BFGS-B wrapper + multi-start
src/plot/       All matplotlib visualization
src/metrics/    Post-run scoring
scripts/        Entry points (run_optimization.py)
data/element_patterns/   CST exports, one .txt per element
results/        Timestamped output folders
tests/          pytest checks
```

Module boundaries are strict: no plotting in `optimize/`, no optimization in `plot/`, etc.

## Hard Rules

1. **Ask before assuming.** If a design decision is ambiguous (angle convention,
   default polarization, missing config key, file format variant), list the
   options with trade-offs and wait for confirmation. Do not pick silently.

2. **Describe before coding.** For any new feature, state the intended approach
   in 2–4 sentences and wait for approval.

3. **Parser first.** When a new CST file format appears, validate the parser
   output (print shapes, sample values) before touching any other module.

4. **No silent defaults.** Missing required config keys → raise a descriptive
   error naming the key. Never fall back to a hardcoded value.

5. **MATLAB portability.** This code will be ported to MATLAB. Mark any
   Python-specific construct with `# [MATLAB] ...` inline. Details in `docs/STYLE.md`.

6. **Update the session log.** At the end of every working session, append an
   entry to `docs/notes.md` under `## Session Log` (template is in that file).

## Physics Conventions (non-negotiable)

- Internal computation uses **radians**; display and CSV output use **degrees**.
- Optimization operates on **linear field magnitude**; dB is for display only.
- Element 0 is the reference element (no normalization unless requested).
- Default optimization target is **co-polarization** (`Abs(Copol)` + `Phase(Copol)`).
  Switch to total field only when `polarization: "total"` is set in `config.yaml`.

## Tooling

- Python: numpy, scipy.optimize (L-BFGS-B), matplotlib, pyyaml
- Tests: pytest
- No linter/formatter is configured yet — follow `docs/STYLE.md` manually

## Out of Scope (do not implement unless asked)

- Mutual coupling correction
- Dynamic jammer tracking (future project, separate codebase)
