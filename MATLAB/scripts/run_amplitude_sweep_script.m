% RUN_AMPLITUDE_SWEEP_SCRIPT  How performance depends on signal vs jammer amplitude.
%
%   A 2-D parameter study over the two power knobs the anti-jam engine exposes,
%   both in dB relative to its fixed sigma_n^2 = 1 per-element noise floor
%   (sim_engine_init):
%
%     y axis  antijam.sigma_s_db   desired-signal power
%     x axis  scenario jn_ratio_db jammer-to-noise ratio
%
%   Because both are referenced to the SAME noise floor, the jammer-to-signal
%   ratio J/S = jn_ratio_db - sigma_s_db falls out as a family of 45-degree
%   lines on the plane (drawn as contours on every heatmap). That is what makes
%   this the right plane to sweep rather than J/S alone: it separates the two
%   distinct ways the link fails — signal-limited (bottom edge: SINR below
%   threshold because the signal is weak, jammer irrelevant) and jammer-limited
%   (bottom-right: SINR below threshold because the nuller cannot keep up) —
%   which a single J/S axis collapses together.
%
%   Per cell it runs closed_loop_run exactly as run_mode_c_demo_script does and
%   scores:
%       availability [%]           fraction of the run with SINR >= sinr_min_db
%       dead time [s]              total time below that threshold
%       mean / steady-state SINR   [dB]
%       oracle gap [dB]            shortfall vs the perfect-knowledge LCMV
%       directivity toward target  [dBi] at (theta_s, phi_s) — "peak gain to
%                                  the target location", via
%                                  compute_directivity_trace (same normalizer
%                                  as the pattern heatmaps, so the number is
%                                  directly comparable to those figures)
%
%   THREE SWEEP DIMENSIONS, each producing its own figure set:
%     * scenario     static always-on / on-off / drift (SWEEP CONFIG below)
%     * loading mode 'adaptive' (P9 adapt.loading_factor_db, data-driven) vs
%                    'fixed'    (P2 adapt.diagonal_loading_db alone).
%                    This plane is exactly what motivated P9 — the claim is
%                    that one untuned formula covers both power regimes — so
%                    the sweep doubles as its validation. A difference figure
%                    (fixed - adaptive) is produced per scenario.
%     * algorithm    oracle (perfect-knowledge upper bound) + lcmv (P2 reactive
%                    covariance tracker). Add 'predict' below at the cost of
%                    roughly tripling the runtime.
%
%   READING THE MAPS — one caveat worth knowing before you look. Availability
%   is scored against a FIXED sinr_min_db, so raising sigma_s_db raises SINR
%   for free: expect availability to saturate at 100% over most of the upper
%   half and collapse across a diagonal cliff, rather than varying smoothly.
%   Oracle gap and directivity are the metrics that carry real information
%   about adaptation QUALITY across the plane; availability and dead time tell
%   you where the operating point sits relative to the threshold. All are
%   plotted for that reason.
%
%   Outputs land in results/amplitude_sweep/<timestamp>/:
%       sweep_<SCN>_adaptive.png   6 metric heatmaps, adaptive loading
%       sweep_<SCN>_fixed.png      6 metric heatmaps, fixed loading
%       sweep_<SCN>_loading.png    fixed vs adaptive + difference maps
%       sweep_<SCN>_oracle.png     the perfect-knowledge reference plane
%       amplitude_sweep.csv        every cell, tidy/long format
%       amplitude_sweep.mat        full `sweep` struct, for re-plotting without
%                                  re-simulating (see REPLOT at the bottom)
%       sweep_params.txt           what was swept
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6, P9].

clear; clc;

% ── 0. Paths: reuse the Milestone-1 port + the anti-jam library ─────
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(fileparts(script_dir));           % <repo> = MATLAB/../
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

% ══════════════════════════════════════════════════════════════════
%  SWEEP CONFIG — everything you would normally want to edit is here
% ══════════════════════════════════════════════════════════════════

% Amplitude grids [dB re noise floor]. Must be uniformly spaced (imagesc).
% 0:5:30 on both axes spans J/S from -30 to +30 dB, covering both regimes the
% milestone has been tuned against: the P6-calibrated weak-signal point
% (sigma_s ~ 0-5, jn 20 -> J/S ~ +15 dB) and the P9 regression point
% (sigma_s 30, jn 20 -> J/S = -10 dB). Widen or refine freely — cost is
% linear in the number of cells and one cell is ~2 s.
sigma_s_db_grid  = 0:5:30;
jn_ratio_db_grid = 0:5:30;

% Algorithms. 'oracle' is mandatory (it defines the oracle-gap reference) and
% is unaffected by the loading mode — it builds R analytically and calls
% adapt_lcmv with zero loading — so it is run ONCE per cell, not once per mode.
% Append 'predict' for the P8 anticipatory nuller (much slower: MUSIC eig on
% the full far-field grid every step).
algorithms = {'lcmv'};

% Diagonal-loading modes to compare (P9). 'adaptive' needs
% adapt.loading_factor_db in config.yaml; 'fixed' removes that key, which is
% exactly how adapt_tracking_init falls back to adapt.diagonal_loading_db.
loading_modes = {'adaptive', 'fixed'};

% Scenarios. The jammer position is FIXED (not the usual per-seed random draw)
% so that amplitude is the only thing varying across cells — otherwise each
% cell would face a different geometry and the maps would be unreadable.
% (90, 200) sits 60 deg from the configured target at (90, 260), well outside
% the 5 deg guard cap. jn_ratio_db is injected per cell by the sweep loop.
sweep_scenarios = { ...
    struct('id', 'STATIC', 'motion', 'static', 'power', 'constant', ...
           'theta_j_deg', 90.0, 'phi_j_deg', 200.0, 'duration_s', 60.0), ...
    struct('id', 'ONOFF',  'motion', 'static', 'power', 'onoff', ...
           'theta_j_deg', 90.0, 'phi_j_deg', 200.0, ...
           'duty_cycle', 0.5, 'toggle_period_s', 20.0, 'duration_s', 100.0), ...
    struct('id', 'DRIFT',  'motion', 'drift',  'power', 'constant', ...
           'theta_j_deg', 90.0, 'phi_j_deg', 200.0, ...
           'theta_drift_deg_per_s', 2.0, 'phi_drift_deg_per_s', 0.0, ...
           'duration_s', 60.0)};

% Fraction of the run treated as "steady state" for the _ss metrics: the last
% (1 - ss_start_frac) of the timeline, so initial convergence is excluded.
ss_start_frac = 0.5;

% Monte Carlo seeds per cell. 1 = the config seed only. Raising this averages
% out seed noise in availability/dead time at a proportional runtime cost.
n_seeds = 1;

% kpi_evaluate's null-pointing-error KPI scans the full (theta, phi) grid for
% local minima at EVERY step — measured at ~4.5 s per run, which is ~95% of
% the total cost of a cell and produces a metric this sweep does not map. The
% metrics below are therefore computed inline (sweep_metrics, at the bottom,
% mirroring the definitions in kpi_evaluate / plot_mode_c_comparison verbatim).
% Set true to call the authoritative kpi_evaluate instead — same numbers for
% everything plotted here, plus null-pointing error in the .mat, ~20x slower.
full_kpi = false;

% ══════════════════════════════════════════════════════════════════

% ── 1. Config: array, polarization, adapt settings ─────────────────
config = read_config_yaml(fullfile(repo_root, 'config.yaml'));
aj_base = config.antijam;

% Loading-mode variants of the adapt section, built once.
adapt_variants = struct();
for i = 1:numel(loading_modes)
    switch loading_modes{i}
        case 'adaptive'
            if ~isfield(config.adapt, 'loading_factor_db') || isempty(config.adapt.loading_factor_db)
                error('run_amplitude_sweep:NoLoadingFactor', ...
                    ['Loading mode ''adaptive'' needs adapt.loading_factor_db in ' ...
                     'config.yaml (the P9 opt-in key); it is absent or empty. Set ' ...
                     'it, or drop ''adaptive'' from loading_modes.']);
            end
            adapt_variants.adaptive = config.adapt;
        case 'fixed'
            a = config.adapt;
            if isfield(a, 'loading_factor_db')
                a = rmfield(a, 'loading_factor_db');   % -> adapt_tracking_init falls
            end                                        %    back to diagonal_loading_db
            adapt_variants.fixed = a;
        otherwise
            error('run_amplitude_sweep:BadLoadingMode', ...
                'Unknown loading mode ''%s'' (expected ''adaptive'' or ''fixed'').', ...
                loading_modes{i});
    end
end

% ── 2. Element patterns -> complex far-field stacks ────────────────
fprintf('Loading element patterns...\n');
patterns  = load_element_patterns(fullfile(repo_root, config.element_patterns_dir));
theta_deg = patterns(1).theta_deg;
phi_deg   = patterns(1).phi_deg;
[stack1, stack2, pol] = select_polarization_stacks(patterns, config);
fprintf('Array: %d elements, grid %dx%d, polarization %s\n', ...
    size(stack1, 1), numel(theta_deg), numel(phi_deg), pol);
fprintf('Target (theta_s=%.0f, phi_s=%.0f) deg; SINR threshold %.1f dB\n', ...
    aj_base.theta_s_deg, aj_base.phi_s_deg, aj_base.sinr_min_db);

% ── 3. Output folder ───────────────────────────────────────────────
timestamp  = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
output_dir = fullfile(repo_root, 'results', 'amplitude_sweep', timestamp);
if ~isfolder(output_dir), mkdir(output_dir); end
fprintf('Output: %s\n', output_dir);

n_sigma = numel(sigma_s_db_grid);
n_jn    = numel(jn_ratio_db_grid);
n_scn   = numel(sweep_scenarios);
n_cells = n_sigma * n_jn * n_scn;
n_runs  = n_cells * n_seeds * (1 + numel(algorithms) * numel(loading_modes));
fprintf('Sweep: %d x %d cells x %d scenarios x %d seed(s) = %d closed-loop runs\n', ...
    n_sigma, n_jn, n_scn, n_seeds, n_runs);

% ── 4. The sweep ───────────────────────────────────────────────────
% sweep.scn{is}.oracle / .alg.<algorithm>.<loading_mode> each hold a struct of
% (n_sigma x n_jn) metric maps, seed-averaged.
sweep = struct();
sweep.sigma_s_db_grid  = sigma_s_db_grid;
sweep.jn_ratio_db_grid = jn_ratio_db_grid;
sweep.algorithms       = algorithms;
sweep.loading_modes    = loading_modes;
sweep.scenario_ids     = cell(1, n_scn);
sweep.scn              = cell(1, n_scn);
sweep.sinr_min_db      = aj_base.sinr_min_db;
sweep.polarization     = pol;
sweep.n_seeds          = n_seeds;
sweep.ss_start_frac    = ss_start_frac;

csv_rows = {};
t_all    = tic;
i_run    = 0;
for is = 1:n_scn
    scn_cfg = sweep_scenarios{is};
    sweep.scenario_ids{is} = scn_cfg.id;
    cell_store = struct('oracle', empty_maps(n_sigma, n_jn), 'alg', struct());
    for ia = 1:numel(algorithms)
        for il = 1:numel(loading_modes)
            cell_store.alg.(algorithms{ia}).(loading_modes{il}) = ...
                empty_maps(n_sigma, n_jn);
        end
    end

    for iy = 1:n_sigma
        for ix = 1:n_jn
            aj = aj_base;
            aj.sigma_s_db = sigma_s_db_grid(iy);
            cfg_cell = scn_cfg;
            cfg_cell.jn_ratio_db = jn_ratio_db_grid(ix);   % per-scenario override

            seed_acc = struct();     % metric -> vector over seeds, per run key
            for iseed = 1:n_seeds
                sim_cfg      = config.sim;
                sim_cfg.seed = config.sim.seed + iseed - 1;
                scn = sim_scenario(cfg_cell, aj, sim_cfg);

                % Oracle first: its SINR timeline is every other run's reference.
                % Loading mode is irrelevant to it (it builds R analytically and
                % calls adapt_lcmv with zero loading), so it runs once per cell.
                o_log = closed_loop_run('oracle', stack1, stack2, theta_deg, ...
                    phi_deg, scn, aj, sim_cfg, config, []);
                o_log.oracle_sinr_db = o_log.sinr_db;
                i_run = i_run + 1;
                seed_acc = accumulate(seed_acc, 'oracle', sweep_metrics( ...
                    o_log, scn, aj, stack1, stack2, theta_deg, phi_deg, ...
                    ss_start_frac, full_kpi));

                for ia = 1:numel(algorithms)
                    for il = 1:numel(loading_modes)
                        cfg_run = config;
                        cfg_run.adapt = adapt_variants.(loading_modes{il});
                        log = closed_loop_run(algorithms{ia}, stack1, stack2, ...
                            theta_deg, phi_deg, scn, aj, sim_cfg, cfg_run, []);
                        log.oracle_sinr_db = o_log.sinr_db;
                        i_run = i_run + 1;
                        seed_acc = accumulate(seed_acc, ...
                            [algorithms{ia} '__' loading_modes{il}], ...
                            sweep_metrics(log, scn, aj, stack1, stack2, ...
                                theta_deg, phi_deg, ss_start_frac, full_kpi));
                    end
                end
            end

            % Seed-average and file into the metric maps.
            cell_store.oracle = store_cell(cell_store.oracle, iy, ix, ...
                mean_over_seeds(seed_acc.oracle));
            csv_rows{end + 1} = csv_row(scn_cfg.id, 'oracle', '-', ...
                sigma_s_db_grid(iy), jn_ratio_db_grid(ix), ...
                mean_over_seeds(seed_acc.oracle)); %#ok<SAGROW>
            for ia = 1:numel(algorithms)
                for il = 1:numel(loading_modes)
                    key = [algorithms{ia} '__' loading_modes{il}];
                    m   = mean_over_seeds(seed_acc.(key));
                    cell_store.alg.(algorithms{ia}).(loading_modes{il}) = ...
                        store_cell(cell_store.alg.(algorithms{ia}).(loading_modes{il}), ...
                            iy, ix, m);
                    csv_rows{end + 1} = csv_row(scn_cfg.id, algorithms{ia}, ...
                        loading_modes{il}, sigma_s_db_grid(iy), ...
                        jn_ratio_db_grid(ix), m); %#ok<SAGROW>
                end
            end

            elapsed = toc(t_all);
            fprintf('  [%s] sigma_s %+5.1f dB, jn %+5.1f dB  (J/S %+5.1f)  |  %d/%d runs, %.0f s elapsed, ~%.0f s left\n', ...
                scn_cfg.id, sigma_s_db_grid(iy), jn_ratio_db_grid(ix), ...
                jn_ratio_db_grid(ix) - sigma_s_db_grid(iy), i_run, n_runs, ...
                elapsed, elapsed * (n_runs - i_run) / max(i_run, 1));
        end
    end
    sweep.scn{is} = cell_store;
end
fprintf('Sweep complete in %.0f s.\n', toc(t_all));

% ── 5. Persist raw results (CSV + .mat) before plotting ────────────
fid = fopen(fullfile(output_dir, 'amplitude_sweep.csv'), 'w');
fprintf(fid, ['scenario,algorithm,loading_mode,sigma_s_db,jn_ratio_db,js_db,' ...
              'availability_pct,dead_time_s,sinr_mean_db,sinr_ss_db,' ...
              'oracle_gap_mean_db,oracle_gap_ss_db,dir_s_dbi_mean,dir_s_dbi_ss,' ...
              'recovery_mean_steps\n']);
fprintf(fid, '%s\n', csv_rows{:});
fclose(fid);
save(fullfile(output_dir, 'amplitude_sweep.mat'), 'sweep');

fid = fopen(fullfile(output_dir, 'sweep_params.txt'), 'w');
fprintf(fid, 'Amplitude sweep — %s\n\n', timestamp);
fprintf(fid, 'array            : %s (%d elements, polarization %s)\n', ...
    config.element_patterns_dir, size(stack1, 1), pol);
fprintf(fid, 'target           : theta_s %.1f, phi_s %.1f deg\n', ...
    aj_base.theta_s_deg, aj_base.phi_s_deg);
fprintf(fid, 'sinr_min_db      : %.1f\n', aj_base.sinr_min_db);
fprintf(fid, 'sigma_s_db grid  : %s\n', mat2str(sigma_s_db_grid));
fprintf(fid, 'jn_ratio_db grid : %s\n', mat2str(jn_ratio_db_grid));
fprintf(fid, 'algorithms       : %s (+ oracle)\n', strjoin(algorithms, ', '));
fprintf(fid, 'loading modes    : %s\n', strjoin(loading_modes, ', '));
fprintf(fid, 'scenarios        : %s\n', strjoin(cellfun(@(s) s.id, ...
    sweep_scenarios, 'UniformOutput', false), ', '));
fprintf(fid, 'seeds per cell   : %d (base %d)\n', n_seeds, config.sim.seed);
fprintf(fid, 'steady state     : last %.0f%% of each run\n', 100 * (1 - ss_start_frac));
fprintf(fid, 'snapshots/step K : %d\n', config.sim.snapshots_per_step);
fprintf(fid, 'forgetting lambda: %.3f\n', config.adapt.forgetting_lambda);
fprintf(fid, 'total runs       : %d in %.0f s\n', n_runs, toc(t_all));
fclose(fid);

% ── 6. Figures ─────────────────────────────────────────────────────
fprintf('Rendering heatmaps...\n');
for is = 1:n_scn
    sid = sweep.scenario_ids{is};
    cs  = sweep.scn{is};

    % Perfect-knowledge reference plane: what the array can do at all here.
    plot_amplitude_heatmaps(sigma_s_db_grid, jn_ratio_db_grid, ...
        metric_panels(cs.oracle, false), ...
        sprintf('%s — oracle (perfect-knowledge upper bound), %s', sid, pol), ...
        fullfile(output_dir, sprintf('sweep_%s_oracle.png', sid)));

    for ia = 1:numel(algorithms)
        alg = algorithms{ia};
        for il = 1:numel(loading_modes)
            lm = loading_modes{il};
            plot_amplitude_heatmaps(sigma_s_db_grid, jn_ratio_db_grid, ...
                metric_panels(cs.alg.(alg).(lm), true), ...
                sprintf('%s — %s, %s diagonal loading, %s', sid, alg, lm, pol), ...
                fullfile(output_dir, sprintf('sweep_%s_%s_%s.png', sid, alg, lm)));
        end

        % P9 comparison: does data-driven loading actually cover both regimes?
        if all(ismember({'adaptive', 'fixed'}, loading_modes))
            a = cs.alg.(alg).adaptive;
            f = cs.alg.(alg).fixed;
            panels = struct( ...
                'map', {f.oracle_gap_ss_db, a.oracle_gap_ss_db, ...
                        f.oracle_gap_ss_db - a.oracle_gap_ss_db, ...
                        f.availability_pct, a.availability_pct, ...
                        a.availability_pct - f.availability_pct}, ...
                'title', {'oracle gap, FIXED loading', 'oracle gap, ADAPTIVE loading', ...
                          'gap improvement (fixed - adaptive)', ...
                          'availability, FIXED loading', 'availability, ADAPTIVE loading', ...
                          'availability gain (adaptive - fixed)'}, ...
                'cbar_label', {'dB', 'dB', 'dB (>0: adaptive wins)', ...
                               '%', '%', '% (>0: adaptive wins)'}, ...
                'style', {'sequential', 'sequential', 'diverging', ...
                          'sequential', 'sequential', 'diverging'}, ...
                'clim', {[], [], [], [0 100], [0 100], []});
            % Common color scale on the two gap panels so they are comparable.
            g = [panels(1).map(:); panels(2).map(:)];
            g = g(isfinite(g));
            if ~isempty(g)
                panels(1).clim = [min(g), max(g)];
                panels(2).clim = panels(1).clim;
            end
            plot_amplitude_heatmaps(sigma_s_db_grid, jn_ratio_db_grid, panels, ...
                sprintf('%s — %s: P9 data-driven vs fixed diagonal loading, %s', ...
                    sid, alg, pol), ...
                fullfile(output_dir, sprintf('sweep_%s_%s_loading.png', sid, alg)));
        end
    end
end

fprintf('\nDone. All artifacts in:\n  %s\n', output_dir);
fprintf(['\nREPLOT without re-simulating:\n' ...
         '  load(''%s'');\n' ...
         '  plot_amplitude_heatmaps(sweep.sigma_s_db_grid, sweep.jn_ratio_db_grid, ...\n' ...
         '      <panels>, ''title'', ''out.png'');\n'], ...
    fullfile(output_dir, 'amplitude_sweep.mat'));


% ────────────────────────── HELPERS ───────────────────────────────

function m = sweep_metrics(log, scn, aj, stack1, stack2, theta_deg, phi_deg, ...
                           ss_start_frac, full_kpi)
% Scalar performance metrics for one closed-loop run.
%
% Definitions mirror the authoritative ones verbatim:
%   availability / dead time  -> plot_mode_c_comparison
%   oracle gap / recovery     -> kpi_evaluate
%   directivity toward target -> compute_directivity_trace (library function)
% They are recomputed here rather than taken from kpi_evaluate because that
% function's null-pointing KPI scans the whole far-field grid for local minima
% at every step (~95% of a cell's cost) and this sweep does not map it. Set
% full_kpi = true in the SWEEP CONFIG to route through kpi_evaluate instead.
thr  = aj.sinr_min_db;
sinr = log.sinr_db;
t    = scn.t_s;
dt   = t(2) - t(1);
ss   = t >= ss_start_frac * t(end);            % steady-state mask

if full_kpi
    k = kpi_evaluate(log, scn, aj);
    m.availability_pct    = 100 * k.availability;
    m.oracle_gap_mean_db  = k.oracle_gap_mean_db;
    m.recovery_mean_steps = mean(k.recovery_steps, 'omitnan');
    m.null_err_median_deg = median(k.null_pointing_err_deg( ...
        ~isnan(k.null_pointing_err_deg)));
else
    m.availability_pct    = 100 * mean(sinr >= thr);
    m.oracle_gap_mean_db  = mean(log.oracle_sinr_db - sinr);
    m.recovery_mean_steps = recovery_mean(sinr, scn, thr);
end
m.dead_time_s      = dt * sum(sinr < thr);
% NOTE: means are taken in dB (a geometric mean in linear power), matching
% kpi_evaluate's oracle_gap_mean_db convention.
m.sinr_mean_db     = mean(sinr);
m.sinr_ss_db       = mean(sinr(ss));
m.oracle_gap_ss_db = mean(log.oracle_sinr_db(ss) - sinr(ss));

dir_s = compute_directivity_trace(log, stack1, stack2, theta_deg, phi_deg);
m.dir_s_dbi_mean = mean(dir_s);
m.dir_s_dbi_ss   = mean(dir_s(ss));
end


function r = recovery_mean(sinr, scn, thr)
% Mean steps from each scenario event until SINR re-crosses the threshold
% (kpi_evaluate section 2). NaN when the scenario defines no events.
n_ev = numel(scn.events);
if n_ev == 0
    r = NaN;
    return
end
rec = NaN(1, n_ev);
for i = 1:n_ev
    k0 = find(scn.t_s >= scn.events{i}.t_s, 1);
    v  = find(sinr(k0:end) >= thr, 1) - 1;
    if ~isempty(v), rec(i) = v; end
end
r = mean(rec, 'omitnan');
end


function maps = empty_maps(n_sigma, n_jn)
% Allocate one NaN (n_sigma x n_jn) map per metric name.
names = metric_names();
maps  = struct();
for i = 1:numel(names)
    maps.(names{i}) = NaN(n_sigma, n_jn);
end
end


function names = metric_names()
names = {'availability_pct', 'dead_time_s', 'sinr_mean_db', 'sinr_ss_db', ...
         'oracle_gap_mean_db', 'oracle_gap_ss_db', 'dir_s_dbi_mean', ...
         'dir_s_dbi_ss', 'recovery_mean_steps'};
end


function acc = accumulate(acc, key, m)
% Append one seed's metric struct to acc.(key) (a struct of vectors).
names = metric_names();
if ~isfield(acc, key)
    acc.(key) = struct();
    for i = 1:numel(names), acc.(key).(names{i}) = []; end
end
for i = 1:numel(names)
    acc.(key).(names{i})(end + 1) = m.(names{i});
end
end


function m = mean_over_seeds(v)
% Collapse a struct of per-seed vectors to a struct of scalars.
names = metric_names();
m = struct();
for i = 1:numel(names)
    m.(names{i}) = mean(v.(names{i}), 'omitnan');
end
end


function maps = store_cell(maps, iy, ix, m)
% File one seed-averaged metric struct into the (iy, ix) cell of every map.
names = metric_names();
for i = 1:numel(names)
    maps.(names{i})(iy, ix) = m.(names{i});
end
end


function row = csv_row(scn_id, alg, loading_mode, sigma_s_db, jn_ratio_db, m)
row = sprintf('%s,%s,%s,%g,%g,%g,%.2f,%.3f,%.2f,%.2f,%.3f,%.3f,%.2f,%.2f,%.2f', ...
    scn_id, alg, loading_mode, sigma_s_db, jn_ratio_db, jn_ratio_db - sigma_s_db, ...
    m.availability_pct, m.dead_time_s, m.sinr_mean_db, m.sinr_ss_db, ...
    m.oracle_gap_mean_db, m.oracle_gap_ss_db, m.dir_s_dbi_mean, m.dir_s_dbi_ss, ...
    m.recovery_mean_steps);
end


function panels = metric_panels(maps, with_gap)
% The standard 6-panel metric set for one algorithm x loading mode.
% with_gap = false for the oracle (its own gap is identically zero).
panels = struct( ...
    'map',        {maps.availability_pct, maps.dead_time_s, maps.sinr_ss_db}, ...
    'title',      {'SINR availability', 'dead time (below threshold)', ...
                   'steady-state SINR'}, ...
    'cbar_label', {'%', 's', 'dB'}, ...
    'style',      {'sequential', 'sequential', 'sequential'}, ...
    'clim',       {[0 100], [], []});
panels(4) = struct('map', maps.dir_s_dbi_ss, ...
    'title', 'directivity toward target (steady state)', ...
    'cbar_label', 'dBi', 'style', 'sequential', 'clim', []);
panels(5) = struct('map', maps.sinr_mean_db, ...
    'title', 'mean SINR (whole run)', ...
    'cbar_label', 'dB', 'style', 'sequential', 'clim', []);
if with_gap
    panels(6) = struct('map', maps.oracle_gap_ss_db, ...
        'title', 'steady-state oracle gap', ...
        'cbar_label', 'dB (lower is better)', 'style', 'sequential', 'clim', []);
else
    panels(6) = struct('map', maps.recovery_mean_steps, ...
        'title', 'mean recovery time after jammer events', ...
        'cbar_label', 'steps', 'style', 'sequential', 'clim', []);
end
end
