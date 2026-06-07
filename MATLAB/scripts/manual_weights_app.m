function manual_weights_app(config_path)
% ══════════════════════════════════════════════════════════════════
% MANUAL_WEIGHTS_APP
% Entry-point script for the interactive manual weight tuner.
%
% MATLAB port of scripts/manual_weights.py (Python / tkinter).
% Launches ManualWeightsTuner, which provides a live 2-D radiation-
% pattern heatmap, per-element amplitude/phase sliders, a directive
% table, and live metrics.
%
%   manual_weights_app                     % uses <repo>/config.yaml
%   manual_weights_app('my_config.yaml')   % explicit config path
%
% Part of: Antenna Array Pattern Optimization Tool.
% ══════════════════════════════════════════════════════════════════

% ── Path setup ────────────────────────────────────────────────────
% Add MATLAB/ to the path so ManualWeightsTuner and all compute
% functions are reachable when this script is run from scripts/.
this_dir   = fileparts(mfilename('fullpath'));   % .../MATLAB/scripts
matlab_dir = fileparts(this_dir);                % .../MATLAB
repo_root  = fileparts(matlab_dir);              % repo root

if ~contains(path, matlab_dir, 'IgnoreCase', true)
    addpath(matlab_dir);
end

% ── Resolve config path ───────────────────────────────────────────
if nargin < 1 || isempty(config_path)
    config_path = fullfile(repo_root, 'config.yaml');
end
% Resolve a bare filename against the repo root so the call works
% regardless of the current working directory.
if ~isfile(config_path) && isfile(fullfile(repo_root, config_path))
    config_path = fullfile(repo_root, config_path);
end
if ~isfile(config_path)
    error('manual_weights_app:NotFound', ...
        'Config file not found: %s', config_path);
end

% ── Re-launch loop ────────────────────────────────────────────────
% "Load Config…" in the UI sets tuner.next_config_path then closes
% the window; main() re-opens a fresh tuner with the new config.
% This mirrors the while-loop in scripts/manual_weights.py:main().
while true
    tuner = ManualWeightsTuner(config_path);
    tuner.run();                         % blocks until window closed
    if isempty(tuner.next_config_path)
        break;
    end
    config_path = tuner.next_config_path;
end
end
