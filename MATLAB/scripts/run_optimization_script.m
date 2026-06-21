% Driver: add the MATLAB/matlab_utils port folder to the path and run the pipelines.
%
% All port functions live flat in MATLAB/matlab_utils/, so a single non-recursive
% addpath of that folder is enough. Do NOT use addpath(genpath(repo_root)) — that
% recursively adds the whole repo (.git/, results/, data/, ...), which is slow,
% pollutes the path, and triggers "shadows a built-in" warnings from unrelated folders.

this_dir   = fileparts(mfilename('fullpath'));   % ...\MATLAB\scripts
matlab_dir = fileparts(this_dir);                % ...\MATLAB
addpath(fullfile(matlab_dir, 'matlab_utils'));

% Config names are resolved against the repo root by the scripts themselves,
% so bare filenames work regardless of the current folder.
run_optimization('../../config.yaml');         % full optimize pipeline -> results/optimizations/<ts>/
% compare_classical('../../test_config.yaml');   % classical vs optimizer -> results/compare_classical/<ts>/
