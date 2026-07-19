% RUN_JAMMER_DEMO_SCRIPT  Entry point for the jammer-scenario GIF demo.
%
%   Thin wrapper: adds matlab_utils + antijam_utils to the path and runs the
%   scenarios defined in the repo-root jammer_config.yaml. See run_jammer_demo.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone.

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

output_dir = run_jammer_demo();
fprintf('Jammer demo complete: %s\n', output_dir);
