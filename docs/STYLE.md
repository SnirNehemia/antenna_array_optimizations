# STYLE.md — Coding Style Reference

> Detailed style rules for code in this project. `CLAUDE.md` requires that
> you follow this file; consult it whenever you write or edit Python source.
> When something is ambiguous, ask — do not guess.

---

## Section Headers

Every logical section within a file gets a header comment in this exact format:

```python
# ────────────────────────── SECTION NAME ──────────────────────────
```

Use a full-width line (aim for ~65 characters total including the section name).
Examples:

```python
# ────────────────────────── IMPORTS ───────────────────────────────

# ────────────────────────── CONSTANTS ─────────────────────────────

# ────────────────────────── HELPER FUNCTIONS ──────────────────────

# ────────────────────────── MAIN FUNCTION ─────────────────────────
```

At the top of every file, include a module-level header block:

```python
# ══════════════════════════════════════════════════════════════════
# MODULE NAME
# Brief one-line description of what this module does.
#
# Part of: Antenna Array Pattern Optimization Tool
# ══════════════════════════════════════════════════════════════════
```

---

## Docstrings — Google Style

Every function and class must have a Google-style docstring.
No exceptions, even for short helper functions.

```python
def compute_array_factor(weights_complex, element_patterns_complex, theta_rad, phi_rad):
    """Compute the coherent array factor over the full angular grid.

    Sums the weighted element patterns to produce the total complex
    radiated field: AF(theta, phi) = sum_n  w_n * E_n(theta, phi).

    Args:
        weights_complex (np.ndarray): Complex weights, shape (N_elements,).
            Units: dimensionless (amplitude in V/V, phase in radians internally).
        element_patterns_complex (np.ndarray): Complex element patterns,
            shape (N_elements, N_theta, N_phi). Units: V/m.
        theta_rad (np.ndarray): Elevation angle grid, shape (N_theta,). Units: radians.
        phi_rad (np.ndarray): Azimuth angle grid, shape (N_phi,). Units: radians.

    Returns:
        np.ndarray: Complex array factor, shape (N_theta, N_phi). Units: V/m.
    """
```

Rules:
- Always state **units** in the arg description (`Units: radians`, `Units: V/m`, `Units: dimensionless`).
- Always state **shape** for numpy arrays.
- If a return value is a numpy array, state its shape and units.
- One blank line between the summary sentence and the Args block.

---

## Variable Naming

Use **descriptive, unit-explicit names**. Never abbreviate unless the abbreviation
is a standard physics symbol (e.g., `N_elements`, not `ne` or `num`).

| Good ✓                      | Bad ✗         | Reason                        |
|-----------------------------|---------------|-------------------------------|
| `theta_rad`                 | `t`, `theta`  | Unit ambiguity                |
| `theta_deg`                 | `ang`         | Unclear meaning               |
| `phi_rad`                   | `p`           | Single-letter names forbidden |
| `weights_complex`           | `w`, `wts`    | Too abbreviated               |
| `element_pattern_complex`   | `ep`, `pat`   | Too abbreviated               |
| `n_elements`                | `N`, `ne`     | Ambiguous                     |
| `cost_value`                | `c`, `J`      | Use `J` only in comments/docs |
| `null_depth_db`             | `nd`          | Unit must be in name          |
| `copol_magnitude_vm`        | `mag`         | Unit in name for physical qty |
| `copol_phase_rad`           | `phase`       | Unit in name                  |
| `target_theta_deg`          | `target`      | Be specific                   |

Loop indices are the only exception: `i`, `n`, `k` are acceptable when the
loop variable is clearly described in the preceding comment.

---

## Physics Comments vs. Implementation Comments

**Physics steps** — always comment. These explain *what* is being computed
and *why*, grounded in the underlying electromagnetic or signal-processing theory.

```python
# Array factor: coherent superposition of weighted element patterns.
# Each element contributes w_n * E_n(theta, phi) to the total field.
array_factor = np.sum(weights_complex[:, np.newaxis, np.newaxis] * element_patterns_complex, axis=0)
```

```python
# Convert phase from degrees to radians for internal computation.
# CST exports phase in degrees; all internal math uses radians.
copol_phase_rad = copol_phase_deg * np.pi / 180.0
```

**Implementation steps** — comment when non-obvious. Skip boilerplate.

```python
# Reshape flat parsed data into (N_theta, N_phi) grid.
# Rows are ordered Theta-major (Theta varies fastest in the CST export).
copol_magnitude_grid = copol_magnitude_flat.reshape(n_theta, n_phi)
```

Do NOT comment obvious lines:
```python
# BAD: open the file
with open(filepath, 'r') as file_handle:
```

---

## MATLAB Porting Flags

When using a Python construct that has no direct MATLAB equivalent, add an
inline comment **on the same line or immediately above**:

```python
weights_real_imag = np.concatenate([weights_complex.real, weights_complex.imag])
# [MATLAB] use [real(weights_complex), imag(weights_complex)]

theta_indices = np.where(np.abs(theta_deg - target_theta_deg) <= half_width_deg)[0]
# [MATLAB] theta_indices = find(abs(theta_deg - target_theta_deg) <= half_width_deg)

element_patterns_complex = np.stack([parse_cst_file(p) for p in pattern_paths])
# [MATLAB] avoid list comprehension; use a for-loop to build cell array or 3D matrix
```

Tag format: always `# [MATLAB]` (uppercase, in brackets) so it's grep-able:
```bash
grep -rn "\[MATLAB\]" src/
```

---

## Code Structure Rules

### One responsibility per function
Each function does exactly one thing. If a function needs a multi-line
comment block to explain its two phases, split it into two functions.

### No magic numbers
Every numeric constant gets a named variable with a comment:

```python
# Minimum null depth target: patterns below this are treated as nulls.
# Units: dB (relative to peak)
MIN_NULL_DEPTH_DB = -30.0
```

### Explicit axis arguments
Always specify `axis=` in numpy reductions. Never rely on default behavior:

```python
# Good: intent is clear
mean_power = np.mean(np.abs(array_factor) ** 2, axis=(0, 1))

# Bad: which axis?
mean_power = np.mean(np.abs(array_factor) ** 2)
```

### No chained operations in physics math
Break complex expressions into named intermediate variables:

```python
# Good: each step is traceable
copol_complex = copol_magnitude_vm * np.exp(1j * copol_phase_rad)
weighted_pattern = weight * copol_complex
array_factor += weighted_pattern

# Bad: opaque
array_factor += w * mag * np.exp(1j * ph * np.pi / 180)
```

---

## Line Length and Formatting

- Max line length: **100 characters**
- Use 4-space indentation (no tabs)
- One blank line between logical blocks within a function
- Two blank lines between top-level functions
- Import order: standard library → third-party (numpy, scipy, matplotlib) → local modules
