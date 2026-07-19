% RUN_ANTIJAM_SCRIPT  Entry point for the anti-jam closed-loop campaign.
%
%   Thin wrapper: adds matlab_utils + antijam_utils to the path and runs the
%   full campaign from the repo-root config.yaml. See run_antijam.m.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone.

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

output_dir = run_antijam();
fprintf('Anti-jam campaign complete: %s\n', output_dir);
