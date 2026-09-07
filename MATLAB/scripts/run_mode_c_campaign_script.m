% ==================================================================
%  run_mode_c_campaign_script.m -- [P12b] Mode C stack evaluation campaign
%
%  WHY THIS EXISTS
%  run_acceptance_grid_script measures ONE configuration of the Mode C stack
%  per invocation, and it is a script, so a caller cannot loop over it without
%  its `clearvars` wiping the loop. This wrapper drives that script once per
%  ARM through its GRID_OVERRIDE hook, from inside a function scope where the
%  clearvars is harmless. Every arm therefore re-uses the grid script's exact
%  case-building, preflight and scoring code rather than a copy that can drift.
%
%  AN ARM is one (algorithm set, adapt config) pairing measured over the whole
%  profile. The arms are defined in ARMS below. Because the grid script seeds
%  every algorithm identically and scores all of them against the same oracle
%  run, a difference between two arms is a difference between the configs and
%  not between the channels they were handed.
%
%  WHAT IT ANSWERS
%    1. Where the whole Mode C stack (lcmv + predict) stands on one metric,
%       across 3 arrays x 5 jammer positions x 3 scenarios.
%    2. The open P9 verdict: data-driven (adaptive) diagonal loading vs the
%       fixed diagonal_loading_db, arbitrated on the P12 oracle-tracking score
%       rather than on nine competing heatmap panels.
%
%  COST WARNING (measured 2026-09-07, not estimated)
%  `predict` runs a MUSIC eigendecomposition over the full far-field grid every
%  step. At doa_stride = 1 that is 42x lcmv on patchs_with_monopoles, 126x on
%  spacing0.6 and 4.6x on ManyDipoles. Striding is cheap in time and expensive
%  in DoA accuracy (spacing0.6 DoA RMSE 1.00 -> 1.45 deg at stride 4 -> 2.30 deg
%  at stride 8), so the campaign runs at stride 1 and pays the time.
%
%  Part of: Antenna Array Pattern Optimization Tool -- anti-jam milestone [P12].
% ==================================================================

function run_mode_c_campaign_script(arm_filter)
% ARM_FILTER : optional cellstr of arm tags to run ({} or absent = all).

if nargin < 1, arm_filter = {}; end

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

% ---- Arm definitions --------------------------------------------
% adapt_override is merged over config.yaml's adapt block for every run of the
% arm; the sentinel '<remove>' DELETES a key, which is how the opt-in adaptive
% loading (adapt.loading_factor_db) is switched off.
%
% NOTE on which is the shipped default: config.yaml currently CARRIES
% loading_factor_db: 0, so 'adaptload' is the as-shipped behaviour and
% 'fixedload' is the counterfactual. P9 never settled which should ship.
%
% The two `cv_` arms differ from their baselines by exactly one thing: the
% adapt.predict.cv block is present, so `predict` steers its null at the
% CV-extrapolated angle while the jammer is moving. config.yaml ships that
% block COMMENTED OUT, so the baseline arms are the as-shipped behaviour and
% the cv arms are the candidate. Building it here rather than reading it from
% config.yaml is deliberate: it removes any dependence on when each arm
% happens to re-read the file.
cv_block = struct('lead_steps', 14, 'q_accel_deg_s2', 0.5, 'r_meas_deg', 1.0, ...
    'gate_deg', 10.0, 'min_track_steps', 20, 'min_speed_deg_s', 0.5, ...
    'max_misses', 20);

repo_root_here = fileparts(fileparts(script_dir));
cfg_here    = read_config_yaml(fullfile(repo_root_here, 'config.yaml'));
predict_cv  = cfg_here.adapt.predict;
predict_cv.cv = cv_block;

ARMS = { ...
    struct('tag', 'adaptload', ...
           'algorithms', {{'lcmv', 'predict'}}, ...
           'adapt_override', struct()), ...
    struct('tag', 'fixedload', ...
           'algorithms', {{'lcmv', 'predict'}}, ...
           'adapt_override', struct('loading_factor_db', '<remove>')), ...
    struct('tag', 'cv_adaptload', ...
           'algorithms', {{'lcmv', 'predict'}}, ...
           'adapt_override', struct('predict', predict_cv)), ...
    struct('tag', 'cv_fixedload', ...
           'algorithms', {{'lcmv', 'predict'}}, ...
           'adapt_override', struct('predict', predict_cv, ...
                                    'loading_factor_db', '<remove>'))};

common = struct();
common.profile         = 'A';
common.scenario_filter = {};        % {} = all of profile A (STATIC/DRIFT/WINDOW)
common.n_seeds         = 3;         % P12 measured 3 as ample (std 0.13-0.18 pp)
common.track_tol_db    = 3.0;
common.pass_pct        = 90.0;

tags = cellfun(@(a) a.tag, ARMS, 'UniformOutput', false);
if ~isempty(arm_filter)
    keep = ismember(tags, arm_filter);
    if ~any(keep)
        error('run_mode_c_campaign:EmptyFilter', ...
            'arm_filter {%s} matches no arm; known arms are {%s}.', ...
            strjoin(arm_filter, ', '), strjoin(tags, ', '));
    end
    ARMS = ARMS(keep);
end

fprintf('Mode C campaign: %d arm(s) -- %s\n\n', numel(ARMS), ...
    strjoin(cellfun(@(a) a.tag, ARMS, 'UniformOutput', false), ', '));

t_all = tic;
for ia = 1:numel(ARMS)
    arm = ARMS{ia};
    fprintf(['\n' repmat('=', 1, 70) '\n  ARM %d/%d: %s\n' repmat('=', 1, 70) '\n'], ...
        ia, numel(ARMS), arm.tag);

    ovr = common;
    ovr.algorithms     = arm.algorithms;
    ovr.adapt_override = arm.adapt_override;
    ovr.out_tag        = arm.tag;

    t_arm = tic;
    run_one_arm(script_dir, ovr);
    fprintf('\n  ARM %s done in %.0f s (campaign %.0f s so far).\n', ...
        arm.tag, toc(t_arm), toc(t_all));
end
fprintf('\nCampaign complete in %.0f s.\n', toc(t_all));
end


% -------------------------------------------------------------------
function run_one_arm(script_dir, ovr)
% Runs the grid script in a THROWAWAY function scope. The grid script begins
% with `clearvars -except grid_override__`, which is why this cannot be inlined
% into the loop above: it would clear the loop counter. Everything this scope
% owns is expendable by the time the script runs.
GRID_OVERRIDE = ovr;                                                  %#ok<NASGU>
run(fullfile(script_dir, 'run_acceptance_grid_script.m'));
end
