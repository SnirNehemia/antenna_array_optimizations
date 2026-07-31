% RUN_MODE_C_DEMO_SCRIPT  Compare the Mode C anti-jam algorithms on one scenario.
%
%   Mode C = the observation mode where the receiver has access to per-element
%   SNAPSHOTS (an N_el x K complex matrix per step) and can therefore estimate
%   the spatial covariance R. All three algorithms below run in Mode C:
%
%     oracle   Perfect-knowledge LCMV/MVDR. Builds R analytically from the TRUE
%              jammer angle (which no real algorithm may see) and computes the
%              optimal weights. This is the UPPER BOUND, drawn for reference.
%              -> adapt_lcmv on the analytic covariance (sim_analytic_covariance).
%
%     lcmv     Reactive covariance tracker (milestone P2). Estimates R from the
%              snapshots with exponential forgetting and recomputes MVDR weights
%              every step. It nulls wherever interference CURRENTLY is, so it
%              lags by ~1/(1-lambda) steps and must re-learn the null every time
%              the jammer turns back on.  -> adapt_tracking_init/_update.
%
%     predict  Predictive / anticipatory nuller (milestone P8, the new method).
%              Runs MUSIC on the same R to get an EXPLICIT jammer angle, tracks
%              the on/off pattern with a periodogram, and PRE-FORMS the null
%              (a hard LCMV constraint at the last-known angle) just before the
%              jammer is predicted to turn back on — so there is no re-learning
%              lag at turn-on.  -> adapt_predict_init/_update (MUSIC + FFT).
%
%   The demo runs one on/off jammer scenario through all three, then produces:
%     * vid_<scn>_<alg>.mp4      one video per method: the 2-D radiation pattern
%                                evolving over time, the jammer dot, and the live
%                                SINR trace (save_run_gif);
%     * mode_c_comparison_<scn>.png  SINR-vs-time for all methods + oracle +
%                                threshold, plus availability / dead-time /
%                                oracle-gap bar charts (plot_mode_c_comparison);
%     * doa_diag_<scn>.png       how 'predict' works: MUSIC angle tracking,
%                                presence detection, and the on/off periodogram
%                                (plot_doa_waterfall);
%     * mode_c_stats.txt         the printed statistics table.
%
%   Everything lands in results/mode_c_demo/<timestamp>/.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

clear; clc;

% ── 0. Paths: reuse the Milestone-1 port + the anti-jam library ─────
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(fileparts(script_dir));           % <repo> = MATLAB/../
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

% Timing: record the wall-clock cost of each major stage. tlog() prints the
% stage time live and appends it to `timings`; a summary is printed at the end.
timings = {};
t_all   = tic;

% ── 1. Config: array, polarization, and the adapt/predict settings ─
% We reuse config.yaml so the demo tracks the real-array configuration and the
% same forgetting factor / diagonal loading / predict keys used elsewhere.
t_stage = tic;
config = read_config_yaml(fullfile(repo_root, 'config.yaml'));
aj     = config.antijam;                                 % theta_s, sinr_min_db, guard, ...
timings = tlog('load config', t_stage, timings);

% ── 2. Element patterns -> complex far-field stacks ─────────────────
t_stage   = tic;
patterns  = load_element_patterns(fullfile(repo_root, config.element_patterns_dir));
theta_deg = patterns(1).theta_deg;                       % (1 x N_theta) elevation grid
phi_deg   = patterns(1).phi_deg;                         % (1 x N_phi)   azimuth grid
[stack1, stack2, pol] = select_polarization_stacks(patterns, config);
timings = tlog('load patterns + stacks', t_stage, timings);
fprintf('Array: %d elements, grid %dx%d, polarization %s\n', ...
    size(stack1, 1), numel(theta_deg), numel(phi_deg), pol);
fprintf('Steer to (theta_s=%.0f, phi_s=%.0f) deg; SINR threshold %.0f dB\n', ...
    aj.theta_s_deg, aj.phi_s_deg, aj.sinr_min_db);

% ── 3. Scenario: an ON/OFF jammer at a fixed off-beam direction ─────
% On/off with a 5 s off-gap is the case that separates the methods: lcmv's
% covariance memory (~2.5 s at lambda = 0.98) has fully decayed by the next
% turn-on, so it re-learns the null each time; 'predict' learns the 10 s toggle
% period and pre-nulls. Edit these fields to try other jammer behaviours.
scn_cfg = struct( ...
    'id',              'ONOFF', ...
    'motion',          'static', ...     % fixed angle (see plan S5/S6 for drift)
    'theta_j_deg',     25.0, ...         % jammer elevation [deg]
    'phi_j_deg',       150.0, ...        % jammer azimuth [deg] (opposite the beam)
    'power',           'onoff', ...      % periodic on/off
    'duty_cycle',      0.5, ...           % 50% on
    'toggle_period_s', 10.0, ...         % 5 s ON / 5 s OFF
    'jn_ratio_db',     1.0, ...         % JAMMER power vs the noise floor [dB].
    ...                                   %   Desired-signal power is
    ...                                   %   antijam.sigma_s_db (config.yaml), so
    ...                                   %   jammer-to-signal J/S = jn_ratio_db -
    ...                                   %   sigma_s_db. Raise this (or lower
    ...                                   %   sigma_s_db) to make the jammer stronger.
    'duration_s',      40.0);           % 10 toggle cycles
t_stage = tic;
scn = sim_scenario(scn_cfg, aj, config.sim);
timings = tlog('build scenario', t_stage, timings);
fprintf('Scenario %s: jammer (theta=%.0f, phi=%.0f), toggle %.0f s, %d steps\n', ...
    scn.id, scn.theta_j_deg(1), scn.phi_j_deg(1), scn_cfg.toggle_period_s, numel(scn.t_s));

% ── 4. Output folder (timestamped) + video settings ────────────────
timestamp  = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
output_dir = fullfile(repo_root, 'results', 'mode_c_demo', timestamp);
if ~isfolder(output_dir), mkdir(output_dir); end
fprintf('Output: %s\n', output_dir);
gif_cfg = struct('format', 'mp4', 'max_frames', 150, 'fps', 12, 'dynamic_range_db', 40);

% ── 5. Oracle reference run (perfect-knowledge upper bound) ─────────
% Run first so every other method can be scored against its SINR timeline.
fprintf('\nRunning oracle (reference)...\n');
t_stage = tic;
o_log = closed_loop_run('oracle', stack1, stack2, theta_deg, phi_deg, ...
    scn, aj, config.sim, config, []);
o_log.oracle_sinr_db = o_log.sinr_db;                    % (its own SINR is the oracle)
timings = tlog('closed-loop oracle', t_stage, timings);

% ── 6. Run each Mode C method, score it, and render its video ───────
methods = {'oracle', 'lcmv', 'predict'};                 % order shown in plots
runs    = {};
for m = 1:numel(methods)
    alg = methods{m};
    fprintf('Running %s...\n', alg);
    if strcmp(alg, 'oracle')
        log = o_log;                                     % already computed
    else
        t_stage = tic;
        log = closed_loop_run(alg, stack1, stack2, theta_deg, phi_deg, ...
            scn, aj, config.sim, config, []);
        log.oracle_sinr_db = o_log.sinr_db;              % attach the reference
        timings = tlog(['closed-loop ' alg], t_stage, timings);
    end

    % KPIs (availability, recovery, oracle gap, ...) from the run log.
    kpi = kpi_evaluate(log, scn, aj);
    runs{end + 1} = struct('algorithm', alg, 'run_log', log, ...
        'kpi', kpi, 'scenario', scn); %#ok<SAGROW>

    % Video: 2-D pattern + jammer dot + live SINR trace for this method.
    vid_path = fullfile(output_dir, sprintf('vid_%s_%s.%s', scn.id, alg, gif_cfg.format));
    t_stage = tic;
    save_run_gif(log, scn, stack1, stack2, theta_deg, phi_deg, aj, gif_cfg, ...
        sprintf('%s / %s', scn.id, alg), vid_path);
    timings = tlog(['render video ' alg], t_stage, timings);
end

% ── 7. Comparison figure + statistics table ────────────────────────
cmp_path = fullfile(output_dir, sprintf('mode_c_comparison_%s.png', scn.id));
t_stage = tic;
stats = plot_mode_c_comparison(runs, aj, cmp_path);
timings = tlog('comparison figure', t_stage, timings);

% Print (and save) the statistics table — the "dead time and such".
lines = {};
lines{end+1} = sprintf('Mode C comparison — scenario %s', scn.id);
lines{end+1} = sprintf('%-9s %10s %12s %12s %14s', ...
    'algorithm', 'avail[%]', 'deadtime[s]', 'oraclegap[dB]', 'recovery[steps]');
for i = 1:numel(stats)
    s = stats(i);
    lines{end+1} = sprintf('%-9s %10.1f %12.2f %12s %14s', s.algorithm, ...
        s.availability, s.dead_time_s, ...
        num_or_dash(s.oracle_gap_db), num_or_dash(s.recovery_mean_steps)); %#ok<SAGROW>
end
report = strjoin(lines, newline);
fprintf('\n%s\n', report);
fid = fopen(fullfile(output_dir, 'mode_c_stats.txt'), 'w');
fprintf(fid, '%s\n', report);
fclose(fid);

% ── 8. How 'predict' works: MUSIC + periodogram diagnostics ────────
predict_log = runs{strcmp(methods, 'predict')}.run_log;
doa_path = fullfile(output_dir, sprintf('doa_diag_%s.png', scn.id));
t_stage = tic;
plot_doa_waterfall(predict_log, scn, aj, doa_path);
timings = tlog('doa figure', t_stage, timings);

% ── Timing summary ─────────────────────────────────────────────────
fprintf('\nTiming summary (major stages):\n');
for i = 1:size(timings, 1)
    fprintf('  %-22s %8.2f s\n', timings{i, 1}, timings{i, 2});
end
fprintf('  %-22s %8.2f s\n', 'TOTAL', toc(t_all));

fprintf('\nDemo complete. All artifacts in:\n  %s\n', output_dir);


% ────────────────────────── HELPERS ───────────────────────────────

function timings = tlog(label, t_stage, timings)
% Record and print the elapsed time of one stage; append {label, seconds}.
dt = toc(t_stage);
timings(end + 1, :) = {label, dt};
fprintf('  [time] %-22s %8.2f s\n', label, dt);
end


function s = num_or_dash(x)
% Format a scalar, or '-' when it is NaN (e.g. oracle has no oracle-gap).
if isnan(x), s = '-'; else, s = sprintf('%.2f', x); end
end
