function optimizer_app(config_path)
% ══════════════════════════════════════════════════════════════════
% OPTIMIZER_APP
% Entry-point script for the interactive manual weight tuner.
%
%   optimizer_app                          % uses <repo>/matlab_config.yaml
%   optimizer_app('my_config.yaml')        % explicit config path
%
% Part of: Antenna Array Pattern Optimization Tool.
% ══════════════════════════════════════════════════════════════════

% ── Version gate ──────────────────────────────────────────────────
% This app is validated on MATLAB R2020a (9.8) and newer. Older releases
% lack GUI properties/functions used here (e.g. clim, GridLineWidth,
% uilabel WordWrap) and fail with cryptic "Unrecognized property" errors.
if verLessThan('matlab', '9.8')
    error('optimizer_app:UnsupportedMATLAB', ...
        ['This app requires MATLAB R2020a (9.8) or newer. ', ...
         'Detected: %s.'], version);
end

% ── Path setup ────────────────────────────────────────────────────
this_dir   = fileparts(mfilename('fullpath'));   % .../MATLAB/scripts
matlab_dir = fileparts(this_dir);                % .../MATLAB
repo_root  = fileparts(matlab_dir);              % repo root
utils_dir  = fullfile(matlab_dir, 'matlab_utils');

if ~contains(path, utils_dir, 'IgnoreCase', true)
    addpath(utils_dir);
end

% ── Resolve config path ───────────────────────────────────────────
if nargin < 1 || isempty(config_path)
    % Use matlab_config.yaml as the default config filename
    config_path = fullfile(this_dir, 'matlab_config.yaml');
end
% Resolve a bare filename against the repo root so the call works
% regardless of the current working directory.
if ~isfile(config_path) && isfile(fullfile(repo_root, config_path))
    config_path = fullfile(repo_root, config_path);
end
if ~isfile(config_path)
    error('optimizer_app:NotFound', ...
        'Config file not found: %s', config_path);
end

% ── Re-launch loop ────────────────────────────────────────────────
while true
    tuner = ManualWeightsTuner(config_path);
    tuner.run();                         % blocks until window closed
    if isempty(tuner.next_config_path)
        break;
    end
    config_path = tuner.next_config_path;
end
end
