# Offline Setup Guide

Instructions for installing the antenna array optimization tool on a machine
without internet access (air-gapped lab computer, secure network, etc.).

The process has two phases:
1. **Online machine** — download all installers and wheels.
2. **Offline machine** — install from the downloaded files.

---

## Phase 1 — Online Machine (any Windows 11 x64 machine with internet)

### 1.1 Download the Python installer

Download **Python 3.13.2** (64-bit) from the official site:

```
https://www.python.org/downloads/release/python-3132/
```

File to grab: `python-3.13.2-amd64.exe`

Save it to a transfer folder, e.g. `D:\offline_bundle\`.

### 1.2 Download all package wheels

Open a command prompt and run:

```cmd
pip download -r requirements.txt --dest D:\offline_bundle\wheels\
```

This downloads `.whl` (wheel) and `.tar.gz` files for every package pinned in
`requirements.txt` and all their transitive dependencies. Typical download size
is 80–120 MB.

> **Important — match the target platform.**
> The wheels must match the *offline* machine's OS and CPU architecture, not the
> online machine. If both machines are Windows 11 64-bit Python 3.13, the default
> `pip download` works. If the offline machine runs a different OS or Python
> version, add platform flags:
>
> ```cmd
> pip download -r requirements.txt --dest D:\offline_bundle\wheels\ ^
>     --platform win_amd64 ^
>     --python-version 313 ^
>     --implementation cp ^
>     --abi cp313 ^
>     --only-binary :all:
> ```

### 1.3 Copy the project files

Copy the entire project folder to the transfer medium (USB drive, DVD, etc.):

```
offline_bundle/
  python-3.13.2-amd64.exe   ← Python installer
  wheels/                   ← all .whl files
  anetnna_array_optimizations/  ← project source (or a zip of it)
```

Transfer to the offline machine.

---

## Phase 2 — Offline Machine

### 2.1 Install Python

Run `python-3.13.2-amd64.exe` and follow the installer.

Recommended options:
- Check **"Add python.exe to PATH"** before clicking Install Now.
- Install for all users if multiple people will use the machine.

Verify:

```cmd
python --version
```

Expected output: `Python 3.13.2`

### 2.2 Install packages from local wheels

```cmd
pip install --no-index --find-links D:\offline_bundle\wheels\ -r requirements.txt
```

- `--no-index` prevents pip from trying to reach PyPI.
- `--find-links` points pip at the local folder of downloaded wheels.

Verify the install:

```cmd
python -c "import numpy, scipy, matplotlib, yaml, PIL; print('All OK')"
```

### 2.3 Verify tkinter (for `manual_weights.py`)

tkinter ships with the Python installer. Confirm it is present:

```cmd
python -m tkinter
```

A small test window should appear. If it does not, re-run the Python installer
and make sure **"tcl/tk and IDLE"** is checked under Optional Features.

---

## Phase 3 — Run the tool

Copy the project folder to any location on the offline machine, then run from
the project root:

```cmd
cd path\to\anetnna_array_optimizations
python scripts\run_optimization.py
```

All features work fully offline — no network calls are made at runtime.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `pip install` fails with "no matching distribution" | Wheel platform mismatch | Re-download with explicit `--platform` flags (see Phase 1.2) |
| `import tkinter` fails | Tk not installed with Python | Re-run Python installer, enable tcl/tk |
| `ModuleNotFoundError: PIL` | Pillow not installed | Check that `pillow-*.whl` is in the wheels folder and re-run Phase 2.2 |
| `No module named yaml` | PyYAML not installed | Same as above for `PyYAML-*.whl` |
| Pattern GIF not saved | Pillow missing or wrong version | Pillow must be ≥ 9.0 for `save_pattern_gif`; check `pip show pillow` |

---

## Keeping wheels up to date

When `requirements.txt` is updated on an internet-connected machine, re-run
Phase 1.2 to refresh the wheels folder, then repeat Phase 2.2 on the offline
machine.

To check what is currently installed on the offline machine:

```cmd
pip list
pip show numpy scipy matplotlib PyYAML pillow
```
