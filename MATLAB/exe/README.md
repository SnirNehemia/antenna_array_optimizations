# OptimizerApp (standalone executable)

Standalone build of `MATLAB/scripts/optimizer_app.m` (the interactive manual
weight tuner), produced with MATLAB Compiler (`mcc`). Bundles all of
`MATLAB/matlab_utils/`, the default `matlab_config.yaml`, and the default
`data/Dipole/` element-pattern folder, so it runs without a MATLAB license —
only the MATLAB Runtime is needed.

## Running it

Double-click **`Run_OptimizerApp.bat`**, not `OptimizerApp.exe` directly — the
batch file puts the MATLAB Runtime DLLs on `PATH` first. `OptimizerApp.exe`
will fail with `Could not find version 26.1 of the MATLAB Runtime` if launched
on its own without that DLL directory on `PATH`.

- On a machine with full **MATLAB R2026a** installed (this build machine), the
  required DLLs live under `<matlabroot>\runtime\win64` — already set in the
  `.bat` file.
- On a machine with only the **MATLAB Runtime R2026a** installed (no full
  MATLAB), edit `MATLAB_RUNTIME_DIR` in `Run_OptimizerApp.bat` to point at that
  installation's `runtime\win64` folder instead. Get the installer from
  https://www.mathworks.com/products/compiler/matlab-runtime.html (must match
  R2026a exactly).

Once running, use the **Load Config** / **Load Data Folder** buttons in the
GUI to point at a different `config.yaml` or element-pattern folder — the
bundled defaults (`data/Dipole/`) are only a starting point.

## Rebuilding after MATLAB source changes

From the repo root, in MATLAB (with MATLAB Compiler licensed):

```matlab
addpath(genpath('MATLAB/matlab_utils'));
mcc('-m', 'MATLAB/scripts/optimizer_app.m', ...
    '-o', 'OptimizerApp', ...
    '-d', 'MATLAB/exe', ...
    '-a', 'MATLAB/scripts/matlab_config.yaml', ...
    '-a', 'data/Dipole');
```

`mcc` auto-detects all `.m` dependencies (including the `ManualWeightsTuner`
class) via the call graph as long as `matlab_utils/` is on the path during the
build — no need to list each function individually. Requires the Optimization
Toolbox (for the "Run Optimization" button's `fmincon` call).
