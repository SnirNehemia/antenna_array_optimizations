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
%   which a single J/S axis collapses together. Whether that separation is
%   actually needed for a given metric is no longer a matter of belief: the
%   J/S curve figures (plot_js_curves) replot every cell against J/S alone, and
%   a metric whose sigma_s rows collapse onto one curve IS a J/S-only metric.
%
%   [2026-08-31] Those curves have now been run (16x16 grid, 5 seeds) and they
%   REFUTE the justification above, while confirming the choice of plane. The
%   45-degree J/S contours are geometrically true but they do not organize the
%   physics: essentially NOTHING here collapses onto a J/S curve. Directivity
%   loss under fixed loading is nearly FLAT in J/S and stratified purely by
%   sigma_s (-3 dB at sigma_s = 0, -12 to -14 dB at sigma_s = 30, at every
%   jammer level including J/N = 0); oracle gap fans out by sigma_s; and
%   availability is pinned at 100% for every sigma_s >= 8 regardless of J/S,
%   so it is a threshold-vs-signal-power effect, not a J/S effect. The real
%   structure on this plane is HORIZONTAL. Keep the 2-D sweep — but keep it
%   because sigma_s is the dominant variable, not because J/S separates two
%   failure regimes.
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
%       directivity LOSS [dB]      the above minus the quiescent beam's
%                                  directivity toward the target. This is the
%                                  desired-signal-cancellation detector: an
%                                  MPDR beamformer fed snapshots that contain
%                                  the desired signal can win on SINR while
%                                  nulling its own signal, which shows up here
%                                  as a large negative number and nowhere else
%                                  in the metric set.
%       operational status         availability >= avail_floor_pct AND
%                                  directivity loss >= -dir_loss_max_db, as a
%                                  4-level pass/fail code. The point of a
%                                  combined mask is that neither metric alone
%                                  catches the "SINR looks fine, the beam is
%                                  destroyed" corner.
%
%   THREE SWEEP DIMENSIONS, each producing its own figure set:
%     * scenario     static always-on / on-off / drift / fast on-off / single
%                    on-window (SWEEP CONFIG below)
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
%   Oracle gap, directivity loss and the operational mask are the metrics that
%   carry real information about adaptation QUALITY across the plane;
%   availability and dead time tell you where the operating point sits relative
%   to the threshold. All are plotted for that reason.
%
%   Outputs land in results/amplitude_sweep/<timestamp>/:
%       scenarios_overview.png       all scenarios' ground-truth jammer
%                                    timelines on one time axis
%       scenario_<SCN>.png           one scenario in detail, events labelled
%       sweep_<SCN>_<alg>_<mode>.png 9 metric heatmaps
%       sweep_<SCN>_oracle.png       the perfect-knowledge reference plane
%       sweep_<SCN>_<alg>_loading.png  fixed vs adaptive + difference maps
%       js_<SCN>.png                 the same metrics as 1-D curves vs J/S
%       trace_<SCN>_<CELL>.png       full time histories for representative
%                                    cells (one per failure regime)
%       amplitude_sweep.csv          every cell, tidy/long format, written
%                                    incrementally as the sweep runs
%       amplitude_sweep.mat          full `sweep` struct, for re-plotting
%                                    without re-simulating (see REPLOT below)
%       sweep_params.txt             what was swept
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone
%   [P6, P9, P11].

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
% 0:2:30 spans J/S from -30 to +30 dB, covering both regimes the milestone has
% been tuned against: the P6-calibrated weak-signal point (sigma_s ~ 0-5,
% jn 20 -> J/S ~ +15 dB) and the P9 regression point (sigma_s 30, jn 20 ->
% J/S = -10 dB). [2026-08-30] Refined from the first sweep's 5 dB steps: the
% availability cliff there was 1-2 cells wide, i.e. the most interesting
% feature on the plane was also the most under-sampled one. Cost is linear in
% the number of cells and one cell-seed-run is ~0.5 s.
% [P12] Pruned from 0:2:30 (16x16 = 256 cells) to 0:5:30 (7x7 = 49). The
% 2026-08-30 campaign established that the plane's structure is essentially
% HORIZONTAL — sigma_s dominates, and the script's own plot_js_curves figures
% refuted the J/S-collapse hypothesis its header used to assert — so a 16x16
% grid was re-measuring a known trend at 5x the cost. The freed budget goes to
% run_acceptance_grid_script, which varies the array and the jammer angle: two
% axes this sweep holds fixed and never tested at all.
%
% Set profile = 'fine' to restore the 16x16 / 5-scenario / 5-seed form. P11's
% reason for refining the grid was real (the availability cliff is 1-2 cells
% wide at 5 dB steps), it is just not what most runs need.
profile = 'coarse';          % 'coarse' | 'fine'

switch profile
    case 'coarse'
        sigma_s_db_grid  = 0:5:30;
        jn_ratio_db_grid = 0:5:30;
    case 'fine'
        sigma_s_db_grid  = 0:2:30;
        jn_ratio_db_grid = 0:2:30;
    otherwise
        error('run_amplitude_sweep:BadProfile', ...
            'profile must be ''coarse'' or ''fine''; got ''%s''.', profile);
end

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
%
% FASTONOFF is new in this run and is here because of a specific finding from
% 2026-08-03: in ONOFF (20 s toggle) at sigma_s_db >= 25 the fixed-loading
% beamformer's directivity toward the target went NEGATIVE (-3.3 to -7.8 dBi)
% while its SINR still read 30+ dB — classic MPDR desired-signal cancellation,
% and it appeared in ONOFF only, not in STATIC or DRIFT. That points at the
% covariance transient around a power toggle rather than at the power level
% itself. FASTONOFF holds everything else equal and raises the toggle rate 4x
% (5 s period instead of 20 s): if the effect is transient-driven it should get
% worse, and if it is a steady-state property of the on-phase it should not
% move. It also makes the recovery-time metric informative — ONOFF gives ~5
% turn-on events per run, FASTONOFF gives ~20.
%
% WINDOW is the other half of the same question. FASTONOFF asks what happens
% when the jammer toggles FASTER; WINDOW asks what happens ONCE, slowly, with
% enough quiet time on both sides to watch the whole lifecycle: 60 s silent ->
% 60 s jamming -> 60 s silent again. It is the only scenario here that emits a
% 'turn_off' event, so it is the only one that can show whether the beam
% RECOVERS — whether a covariance tracker that dug a null (and, at high
% sigma_s, wrecked its own main lobe doing so) climbs back to the quiescent
% beam once the threat stops, or stays deformed. Steady-state scalars cannot
% answer that; this scenario plus its trace figures can.
% [P12] Trimmed to STATIC / DRIFT / WINDOW. FASTONOFF existed to test one
% hypothesis — that desired-signal cancellation was driven by the covariance
% transient around a power toggle — and the 2026-08-30 campaign killed it: the
% effect is purely sigma_s-driven and STATIC is in fact the worst case. ONOFF's
% multi-cycle averaging is better bought with seeds; what a multi-cycle run can
% uniquely answer (does the tracker accumulate state across cycles?) needs the
% recovery reported PER CYCLE, which is a profile-B case in
% run_acceptance_grid_script, not an averaged cell here. WINDOW remains as the
% clean single turn-on / turn-off lifecycle measurement.
sweep_scenarios = { ...
    struct('id', 'STATIC', 'motion', 'static', 'power', 'constant', ...
           'theta_j_deg', 90.0, 'phi_j_deg', 200.0, 'duration_s', 60.0), ...
    struct('id', 'DRIFT',  'motion', 'drift',  'power', 'constant', ...
           'theta_j_deg', 90.0, 'phi_j_deg', 200.0, ...
           'theta_drift_deg_per_s', 2.0, 'phi_drift_deg_per_s', 0.0, ...
           'duration_s', 60.0), ...
    struct('id', 'WINDOW', 'motion', 'static', 'power', 'window', ...
           'theta_j_deg', 90.0, 'phi_j_deg', 200.0, ...
           'on_time_s', 60.0, 'off_time_s', 120.0, 'duration_s', 180.0)};

% Fraction of the run treated as "steady state" for the _ss metrics: the last
% (1 - ss_start_frac) of the timeline, so initial convergence is excluded.
ss_start_frac = 0.5;

% Monte Carlo seeds per cell. [P12] Reduced from 5 to 3, and this time the
% number is measured rather than guessed: the acceptance grid's per-seed
% standard deviations (now emitted here too, as the <metric>_std columns) put
% the median seed-to-seed spread of the headline score at 0.13-0.18 percentage
% points, with a p90 of 1.4-4.2 pp. Against score differences of interest that
% run 20% vs 98%, three seeds is comfortably enough. Re-read the _std columns
% after any change that moves the noise floor.
n_seeds = 3;

% Operational pass/fail mask thresholds (the 'operational status' panel).
% 90% is the availability floor the milestone gates against elsewhere; 3 dB of
% directivity loss is the point past which the beam is no longer meaningfully
% pointed at the target even if the SINR arithmetic still works out.
avail_floor_pct   = 90.0;
dir_loss_max_db   = 3.0;

% [P12] Oracle-tracking tolerance for the headline track_score_pct metric: the
% run is "keeping up" while it is within this many dB of the perfect-knowledge
% LCMV. See kpi_sweep_metrics.
track_tol_db = 3.0;

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
n_el = size(stack1, 1);
fprintf('Array: %d elements, grid %dx%d, polarization %s\n', ...
    n_el, numel(theta_deg), numel(phi_deg), pol);
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

% ── 4. Scenario overview figures (ground truth, no simulation) ──────
% Drawn FIRST, before anything can fail, so that even an aborted sweep leaves
% behind a statement of what it was going to measure. These consume only
% sim_scenario output — the same arrays fed to sim_engine_step and the oracle,
% and never to any adapt_/agent_ algorithm.
fprintf('Rendering scenario overview figures...\n');
overview_scns = cell(1, n_scn);
for is = 1:n_scn
    cfg_ov = sweep_scenarios{is};
    cfg_ov.jn_ratio_db = aj_base.jn_ratio_db;   % nominal level, for the shape
    overview_scns{is} = sim_scenario(cfg_ov, aj_base, config.sim);
end
try
    plot_scenario_overview(overview_scns, aj_base, sprintf( ...
        'Jammer scenarios (ground truth) — shapes at the nominal J/N = %.0f dB', ...
        aj_base.jn_ratio_db), ...
        fullfile(output_dir, 'scenarios_overview.png'));
    for is = 1:n_scn
        plot_scenario_overview(overview_scns(is), aj_base, sprintf( ...
            'Scenario %s — jammer ground truth (nominal J/N = %.0f dB)', ...
            overview_scns{is}.id, aj_base.jn_ratio_db), ...
            fullfile(output_dir, sprintf('scenario_%s.png', overview_scns{is}.id)));
    end
catch err
    warning('run_amplitude_sweep:OverviewFailed', ...
        'Scenario overview figures failed (%s); continuing with the sweep.', ...
        err.message);
end

% ── 5. The sweep ───────────────────────────────────────────────────
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
sweep.avail_floor_pct  = avail_floor_pct;
sweep.dir_loss_max_db  = dir_loss_max_db;
sweep.track_tol_db     = track_tol_db;
sweep.profile          = profile;

% CSV is opened now and written row-by-row rather than at the end: this run is
% hours long, and a crash in hour three should not cost the first two.
csv_path = fullfile(output_dir, 'amplitude_sweep.csv');
csv_fid  = fopen(csv_path, 'w');
if csv_fid < 0
    error('run_amplitude_sweep:CsvOpen', 'Cannot open %s for writing.', csv_path);
end
% [P12] Header and rows are now generated from metric_names() rather than
% hand-listed, so adding a metric cannot silently desynchronise them; each
% metric contributes a mean column and a <metric>_std column.
fprintf(csv_fid, '%s\n', csv_header());

dir_ref_dbi = NaN;      % quiescent-beam directivity toward the target; filled
                        % on the first oracle run (needs the engine's e_s/grid)
n_failed = 0;
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
                % One cell-seed is wrapped as a unit: a single bad cell must not
                % cost an overnight sweep, but it must not pass silently either.
                % A failure leaves NaNs in that cell and prints loudly.
                try
                    scn = sim_scenario(cfg_cell, aj, sim_cfg);

                    % Oracle first: its SINR timeline is every other run's
                    % reference. Loading mode is irrelevant to it (analytic R,
                    % zero loading), so it runs once per cell-seed.
                    o_log = closed_loop_run('oracle', stack1, stack2, theta_deg, ...
                        phi_deg, scn, aj, sim_cfg, config, []);
                    o_log.oracle_sinr_db = o_log.sinr_db;
                    i_run = i_run + 1;

                    if ~isfinite(dir_ref_dbi)
                        dir_ref_dbi = kpi_quiescent_directivity(o_log, n_el, ...
                            stack1, stack2, theta_deg, phi_deg);
                        fprintf('Quiescent-beam directivity toward the target: %.2f dBi (loss reference)\n', ...
                            dir_ref_dbi);
                    end

                    seed_acc = accumulate(seed_acc, 'oracle', kpi_sweep_metrics( ...
                        o_log, scn, aj, stack1, stack2, theta_deg, phi_deg, ...
                        ss_start_frac, full_kpi, dir_ref_dbi, track_tol_db));

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
                                kpi_sweep_metrics(log, scn, aj, stack1, stack2, ...
                                    theta_deg, phi_deg, ss_start_frac, full_kpi, ...
                                    dir_ref_dbi, track_tol_db));
                        end
                    end
                catch err
                    n_failed = n_failed + 1;
                    warning('run_amplitude_sweep:CellFailed', ...
                        ['[%s] sigma_s %g, jn %g, seed %d FAILED: %s (%s). ' ...
                         'That seed is dropped; the cell keeps its other seeds.'], ...
                        scn_cfg.id, sigma_s_db_grid(iy), jn_ratio_db_grid(ix), ...
                        sim_cfg.seed, err.message, err.identifier);
                end
            end

            % Seed-average and file into the metric maps.
            if isfield(seed_acc, 'oracle')
                m_or  = mean_over_seeds(seed_acc.oracle);
                sd_or = std_over_seeds(seed_acc.oracle);
            else
                m_or  = nan_metrics();    % every seed of this cell failed
                sd_or = nan_std();
            end
            cell_store.oracle = store_cell(cell_store.oracle, iy, ix, m_or);
            fprintf(csv_fid, '%s\n', csv_row(scn_cfg.id, 'oracle', '-', ...
                sigma_s_db_grid(iy), jn_ratio_db_grid(ix), m_or, sd_or));
            for ia = 1:numel(algorithms)
                for il = 1:numel(loading_modes)
                    key = [algorithms{ia} '__' loading_modes{il}];
                    if isfield(seed_acc, key)
                        m  = mean_over_seeds(seed_acc.(key));
                        sd = std_over_seeds(seed_acc.(key));
                    else
                        m  = nan_metrics();
                        sd = nan_std();
                    end
                    cell_store.alg.(algorithms{ia}).(loading_modes{il}) = ...
                        store_cell(cell_store.alg.(algorithms{ia}).(loading_modes{il}), ...
                            iy, ix, m);
                    fprintf(csv_fid, '%s\n', csv_row(scn_cfg.id, algorithms{ia}, ...
                        loading_modes{il}, sigma_s_db_grid(iy), ...
                        jn_ratio_db_grid(ix), m, sd));
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
    % Checkpoint after every scenario, for the same reason the CSV streams.
    save(fullfile(output_dir, 'amplitude_sweep.mat'), 'sweep');
end
fclose(csv_fid);
fprintf('Sweep complete in %.0f s (%d failed cell-seeds).\n', toc(t_all), n_failed);

sweep.dir_ref_dbi = dir_ref_dbi;
save(fullfile(output_dir, 'amplitude_sweep.mat'), 'sweep');

fid = fopen(fullfile(output_dir, 'sweep_params.txt'), 'w');
fprintf(fid, 'Amplitude sweep — %s\n\n', timestamp);
fprintf(fid, 'array            : %s (%d elements, polarization %s)\n', ...
    config.element_patterns_dir, n_el, pol);
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
fprintf(fid, 'quiescent dir    : %.2f dBi (directivity-loss reference)\n', dir_ref_dbi);
fprintf(fid, 'operational mask : availability >= %.0f%% AND dir loss <= %.0f dB\n', ...
    avail_floor_pct, dir_loss_max_db);
fprintf(fid, 'failed cell-seeds: %d\n', n_failed);
fprintf(fid, 'total runs       : %d in %.0f s\n', n_runs, toc(t_all));
fclose(fid);

% ── 6. Heatmap figures ─────────────────────────────────────────────
fprintf('Rendering heatmaps...\n');
for is = 1:n_scn
    sid = sweep.scenario_ids{is};
    cs  = sweep.scn{is};

    % Perfect-knowledge reference plane: what the array can do at all here.
    % Same 9 panels as every algorithm figure, deliberately — a fixed layout is
    % what makes the figures comparable by flipping between them.
    safe_plot(@() plot_amplitude_heatmaps(sigma_s_db_grid, jn_ratio_db_grid, ...
        metric_panels(cs.oracle, avail_floor_pct, dir_loss_max_db), ...
        sprintf('%s — oracle (perfect-knowledge upper bound), %s', sid, pol), ...
        fullfile(output_dir, sprintf('sweep_%s_oracle.png', sid))));

    for ia = 1:numel(algorithms)
        alg = algorithms{ia};
        for il = 1:numel(loading_modes)
            lm = loading_modes{il};
            safe_plot(@() plot_amplitude_heatmaps(sigma_s_db_grid, jn_ratio_db_grid, ...
                metric_panels(cs.alg.(alg).(lm), avail_floor_pct, dir_loss_max_db), ...
                sprintf('%s — %s, %s diagonal loading, %s', sid, alg, lm, pol), ...
                fullfile(output_dir, sprintf('sweep_%s_%s_%s.png', sid, alg, lm))));
        end

        % P9 comparison: does data-driven loading actually cover both regimes?
        if all(ismember({'adaptive', 'fixed'}, loading_modes))
            a = cs.alg.(alg).adaptive;
            f = cs.alg.(alg).fixed;
            panels = [ ...
                mk_panel(f.oracle_gap_ss_db,  'oracle gap, FIXED loading', 'dB', 'sequential'), ...
                mk_panel(a.oracle_gap_ss_db,  'oracle gap, ADAPTIVE loading', 'dB', 'sequential'), ...
                mk_panel(f.oracle_gap_ss_db - a.oracle_gap_ss_db, ...
                    'gap improvement (fixed - adaptive)', 'dB (>0: adaptive wins)', 'diverging'), ...
                mk_panel(f.availability_pct,  'availability, FIXED loading', '%', 'sequential', [0 100]), ...
                mk_panel(a.availability_pct,  'availability, ADAPTIVE loading', '%', 'sequential', [0 100]), ...
                mk_panel(a.availability_pct - f.availability_pct, ...
                    'availability gain (adaptive - fixed)', '% (>0: adaptive wins)', 'diverging'), ...
                mk_panel(f.dir_loss_db_ss,    'directivity loss, FIXED loading', 'dB (<0: signal cancelled)', 'diverging'), ...
                mk_panel(a.dir_loss_db_ss,    'directivity loss, ADAPTIVE loading', 'dB (<0: signal cancelled)', 'diverging'), ...
                mk_panel(a.dir_loss_db_ss - f.dir_loss_db_ss, ...
                    'directivity-loss gain (adaptive - fixed)', 'dB (>0: adaptive wins)', 'diverging')];
            % Common color scale on the two gap panels so they are comparable.
            g = [panels(1).map(:); panels(2).map(:)];
            g = g(isfinite(g));
            if ~isempty(g)
                panels(1).clim = [min(g), max(g)];
                panels(2).clim = panels(1).clim;
            end
            safe_plot(@() plot_amplitude_heatmaps(sigma_s_db_grid, jn_ratio_db_grid, panels, ...
                sprintf('%s — %s: P9 data-driven vs fixed diagonal loading, %s', ...
                    sid, alg, pol), ...
                fullfile(output_dir, sprintf('sweep_%s_%s_loading.png', sid, alg))));
        end
    end
end

% ── 7. J/S collapse curves ─────────────────────────────────────────
% Does the plane reduce to one variable? One figure per scenario, every cell
% replotted at x = J/S, colored by its sigma_s row.
fprintf('Rendering J/S curve figures...\n');
for is = 1:n_scn
    sid = sweep.scenario_ids{is};
    cs  = sweep.scn{is};
    series = {};
    labels = {};
    for ia = 1:numel(algorithms)
        for il = 1:numel(loading_modes)
            series{end + 1} = cs.alg.(algorithms{ia}).(loading_modes{il}); %#ok<SAGROW>
            labels{end + 1} = sprintf('%s / %s loading', algorithms{ia}, ...
                loading_modes{il}); %#ok<SAGROW>
        end
    end
    series{end + 1} = cs.oracle;
    labels{end + 1} = 'oracle';

    jp = [ ...
        mk_curve_panel(series, 'availability_pct', 'SINR availability', ...
            '%', avail_floor_pct), ...
        mk_curve_panel(series, 'dead_time_s', 'dead time (below threshold)', ...
            's', []), ...
        mk_curve_panel(series, 'sinr_ss_db', 'steady-state SINR', ...
            'dB', aj_base.sinr_min_db), ...
        mk_curve_panel(series, 'oracle_gap_ss_db', 'steady-state oracle gap', ...
            'dB (lower is better)', 0), ...
        mk_curve_panel(series, 'dir_loss_db_ss', ...
            'directivity loss vs quiescent beam', 'dB (<0: signal cancelled)', 0), ...
        mk_curve_panel(series, 'recovery_mean_steps', ...
            'mean recovery time after jammer events', 'steps', [])];
    safe_plot(@() plot_js_curves(sigma_s_db_grid, jn_ratio_db_grid, jp, labels, ...
        sprintf(['%s — metrics vs J/S alone, %s. Rows that COLLAPSE onto one ' ...
                 'curve are J/S-governed; rows that FAN OUT need the 2-D plane.'], ...
                sid, pol), ...
        fullfile(output_dir, sprintf('js_%s.png', sid))));
end

% ── 8. Representative-cell traces ──────────────────────────────────
% Six scalars per run is what makes the grid readable and is also what hides
% every transient. Re-run a handful of deliberately chosen cells per scenario
% (one per failure regime) at the base seed and dump the full time histories.
fprintf('Rendering representative-cell traces...\n');
for is = 1:n_scn
    scn_cfg = sweep_scenarios{is};
    sid     = scn_cfg.id;
    cs      = sweep.scn{is};
    ref_maps = cs.alg.(algorithms{1}).(loading_modes{end});   % 'fixed' by default
    picks = pick_representative_cells(ref_maps, sigma_s_db_grid, jn_ratio_db_grid);

    for ip = 1:numel(picks)
        iy = picks(ip).iy;
        ix = picks(ip).ix;
        aj = aj_base;
        aj.sigma_s_db = sigma_s_db_grid(iy);
        cfg_cell = scn_cfg;
        cfg_cell.jn_ratio_db = jn_ratio_db_grid(ix);
        sim_cfg = config.sim;
        sim_cfg.seed = config.sim.seed;      % base seed only — this is a picture

        try
            scn   = sim_scenario(cfg_cell, aj, sim_cfg);
            o_log = closed_loop_run('oracle', stack1, stack2, theta_deg, ...
                phi_deg, scn, aj, sim_cfg, config, []);
            traces = struct('label', 'oracle', 'sinr_db', o_log.sinr_db, ...
                'dir_s_dbi', compute_directivity_trace(o_log, stack1, stack2, ...
                    theta_deg, phi_deg), 'is_oracle', true);
            for ia = 1:numel(algorithms)
                for il = 1:numel(loading_modes)
                    cfg_run = config;
                    cfg_run.adapt = adapt_variants.(loading_modes{il});
                    log = closed_loop_run(algorithms{ia}, stack1, stack2, ...
                        theta_deg, phi_deg, scn, aj, sim_cfg, cfg_run, []);
                    traces(end + 1) = struct( ...
                        'label', sprintf('%s / %s', algorithms{ia}, loading_modes{il}), ...
                        'sinr_db', log.sinr_db, ...
                        'dir_s_dbi', compute_directivity_trace(log, stack1, ...
                            stack2, theta_deg, phi_deg), ...
                        'is_oracle', false); %#ok<SAGROW>
                end
            end
            safe_plot(@() plot_cell_traces(scn, aj, traces, dir_ref_dbi, ...
                sprintf(['%s / %s cell — sigma_s = %g dB, J/N = %g dB ' ...
                         '(J/S = %+g dB), seed %d'], sid, picks(ip).name, ...
                        sigma_s_db_grid(iy), jn_ratio_db_grid(ix), ...
                        jn_ratio_db_grid(ix) - sigma_s_db_grid(iy), sim_cfg.seed), ...
                fullfile(output_dir, sprintf('trace_%s_%s.png', sid, picks(ip).name))));
        catch err
            warning('run_amplitude_sweep:TraceFailed', ...
                'Trace figure %s/%s failed: %s', sid, picks(ip).name, err.message);
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

function maps = empty_maps(n_sigma, n_jn)
% Allocate one NaN (n_sigma x n_jn) map per metric name.
names = metric_names();
maps  = struct();
for i = 1:numel(names)
    maps.(names{i}) = NaN(n_sigma, n_jn);
end
end


function names = metric_names()
% [P12] Delegated to the shared list so this sweep and the acceptance grid
% cannot drift apart on what a metric set contains.
names = kpi_sweep_metric_names();
end


function m = nan_metrics()
% An all-NaN metric struct, for a cell whose every seed failed.
names = metric_names();
m = struct();
for i = 1:numel(names)
    m.(names{i}) = NaN;
end
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


function sd = std_over_seeds(v)
% Per-seed standard deviation beside every mean.
%
% [P12] This is the missing half of the seed question. Before it, seeds were
% averaged inline and nothing else was kept, so "is this difference between two
% loading modes real or is it seed noise?" and "how many seeds do we actually
% need?" were both unanswerable from a completed sweep's own output — the only
% way to find out was to run it again. Writing sd beside mu makes both
% answerable from the CSV alone.
names = metric_names();
sd = struct();
for i = 1:numel(names)
    x = v.(names{i});
    x = x(isfinite(x));
    if numel(x) < 2
        % std of a single sample is NaN in MATLAB; report 0 so a one-seed run
        % does not read as a failed measurement.
        sd.(names{i}) = 0;
    else
        sd.(names{i}) = std(x);
    end
end
end


function sd = nan_std()
% All-NaN companion to nan_metrics, for a cell whose every seed failed.
names = metric_names();
sd = struct();
for i = 1:numel(names)
    sd.(names{i}) = NaN;
end
end


function maps = store_cell(maps, iy, ix, m)
% File one seed-averaged metric struct into the (iy, ix) cell of every map.
names = metric_names();
for i = 1:numel(names)
    maps.(names{i})(iy, ix) = m.(names{i});
end
end


function h = csv_header()
names = metric_names();
cols  = {'scenario', 'algorithm', 'loading_mode', 'sigma_s_db', 'jn_ratio_db', 'js_db'};
for i = 1:numel(names)
    cols{end + 1} = names{i};                                         %#ok<AGROW>
    cols{end + 1} = [names{i} '_std'];                                %#ok<AGROW>
end
h = strjoin(cols, ',');
end


function row = csv_row(scn_id, alg, loading_mode, sigma_s_db, jn_ratio_db, m, sd)
names = metric_names();
vals  = cell(1, 2 * numel(names));
for i = 1:numel(names)
    vals{2 * i - 1} = sprintf('%.4g', m.(names{i}));
    vals{2 * i}     = sprintf('%.4g', sd.(names{i}));
end
row = sprintf('%s,%s,%s,%g,%g,%g,%s', ...
    scn_id, alg, loading_mode, sigma_s_db, jn_ratio_db, ...
    jn_ratio_db - sigma_s_db, strjoin(vals, ','));
end


function p = mk_panel(map, ttl, cbar_label, style, clim, clim_floor, cat_labels)
% Build ONE fully-populated panel struct. Every optional field is present (as
% [] when unused) because MATLAB refuses to concatenate structs whose field
% sets differ, and panels are assembled as a struct array.
if nargin < 5, clim       = []; end
if nargin < 6, clim_floor = []; end
if nargin < 7, cat_labels = {}; end
p = struct('map', map, 'title', ttl, 'cbar_label', cbar_label, ...
    'style', style, 'clim', clim, 'clim_floor', clim_floor, ...
    'cat_labels', {cat_labels});
end


function panels = metric_panels(maps, avail_floor_pct, dir_loss_max_db)
% The standard 9-panel metric set. Identical for the oracle and for every
% algorithm/loading variant on purpose: a fixed layout is what lets two figures
% be compared by flipping between them.
%
% clim_floor is what stops a metric that never left its floor from being
% auto-stretched to full scale. The 2026-08-03 sweep's oracle dead-time panel
% read 0.0-0.1 s out of 60 s and, auto-scaled to [0, 0.05], rendered as
% saturated yellow across the whole plane — i.e. it looked like total failure
% while reporting near-perfect behaviour. Same for recovery time, which is
% 0 or 1 steps for the oracle.
panels = [ ...
    mk_panel(maps.track_score_pct, 'oracle-tracking score (headline)', ...
        '% of run within 3 dB of oracle', 'sequential', [0 100]), ...
    mk_panel(maps.availability_pct, 'SINR availability', '%', 'sequential', [0 100]), ...
    mk_panel(maps.dead_time_s, 'dead time (below threshold)', 's', 'sequential', ...
        [], [0 1]), ...
    mk_panel(maps.recovery_mean_steps, 'mean recovery time after jammer events', ...
        'steps', 'sequential', [], [0 5]), ...
    mk_panel(maps.sinr_mean_db, 'mean SINR (whole run)', 'dB', 'sequential'), ...
    mk_panel(maps.sinr_ss_db, 'steady-state SINR', 'dB', 'sequential'), ...
    mk_panel(maps.oracle_gap_ss_db, 'steady-state oracle gap', ...
        'dB (lower is better)', 'sequential', [], [0 1]), ...
    mk_panel(maps.dir_s_dbi_ss, 'directivity toward target (steady state)', ...
        'dBi', 'sequential'), ...
    mk_panel(maps.dir_loss_db_ss, 'directivity loss vs quiescent beam', ...
        'dB (<0: desired signal cancelled)', 'diverging', [], [-1 1]), ...
    mk_panel(operational_code(maps, avail_floor_pct, dir_loss_max_db), ...
        sprintf('operational status (avail >= %.0f%%, dir loss <= %.0f dB)', ...
            avail_floor_pct, dir_loss_max_db), ...
        '', 'categorical', [], [], ...
        {'both fail', 'SINR fail', 'beam fail', 'pass'})];
end


function code = operational_code(maps, avail_floor_pct, dir_loss_max_db)
% 2 * (availability OK) + (directivity OK), so the four codes are
%   0 both fail | 1 SINR fail (beam OK) | 2 beam fail (SINR OK) | 3 pass.
% Code 2 is the whole reason this panel exists: a run that clears the SINR
% threshold while its beam has collapsed passes every other metric here.
ok_a = maps.availability_pct >= avail_floor_pct;
ok_d = maps.dir_loss_db_ss   >= -dir_loss_max_db;
code = 2 * double(ok_a) + double(ok_d);
code(~isfinite(maps.availability_pct) | ~isfinite(maps.dir_loss_db_ss)) = NaN;
end


function p = mk_curve_panel(series, field, ttl, ylab, yref)
% One plot_js_curves panel: the same metric pulled from every series.
p = struct('data', {cellfun(@(s) s.(field), series, 'UniformOutput', false)}, ...
    'title', ttl, 'ylabel', ylab, 'yref', yref);
end


function picks = pick_representative_cells(maps, sigma_s_db_grid, jn_ratio_db_grid)
% Four cells chosen to sit in four different regimes, so the trace figures
% cover the plane's behaviour rather than four samples of the same thing:
%   signal_limited  weakest signal, weakest jammer — failure is the signal, not
%                   the jammer, and the beamformer has nothing to fix.
%   jammer_limited  weakest signal, strongest jammer — the hardest cell.
%   cliff           whichever cell's availability sits closest to 50%, i.e. on
%                   the diagonal transition the heatmaps show as a step. This
%                   is the only one chosen from the DATA rather than from the
%                   grid corners, because its location is what the sweep found.
%   high_snr        strongest signal at a mid jammer level — the corner where
%                   the 2026-08-03 sweep found negative directivity.
n_sigma = numel(sigma_s_db_grid);
n_jn    = numel(jn_ratio_db_grid);
[~, ix_mid] = min(abs(jn_ratio_db_grid - 20));

a = maps.availability_pct;
[~, k_cliff] = min(abs(a(:) - 50));
[iy_c, ix_c] = ind2sub(size(a), k_cliff);
if ~isfinite(a(k_cliff))               % all-NaN map: fall back to the centre
    iy_c = ceil(n_sigma / 2);
    ix_c = ceil(n_jn / 2);
end

picks = struct( ...
    'name', {'signal_limited', 'jammer_limited', 'cliff', 'high_snr'}, ...
    'iy',   {1, 1, iy_c, n_sigma}, ...
    'ix',   {1, n_jn, ix_c, ix_mid});

% Drop duplicates (the cliff can land on a corner) so the same figure is not
% rendered twice under two names.
seen = false(1, numel(picks));
keys = arrayfun(@(p) p.iy * 1000 + p.ix, picks);
for i = 1:numel(picks)
    seen(i) = any(keys(1:i - 1) == keys(i));
end
picks = picks(~seen);
end


function safe_plot(fn)
% Run one figure call, converting a failure into a warning. The data is already
% on disk by the time any of these run, and a graphics-driver hiccup at hour
% three of an overnight sweep must not take the other figures down with it.
try
    fn();
catch err
    warning('run_amplitude_sweep:PlotFailed', ...
        'A figure failed to render (%s): %s', err.identifier, err.message);
end
end
