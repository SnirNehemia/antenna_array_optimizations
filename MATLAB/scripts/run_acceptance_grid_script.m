% ══════════════════════════════════════════════════════════════════
%  run_acceptance_grid_script.m — [P12] tiered anti-jam test suite
%
%  WHY THIS EXISTS
%  The amplitude sweep (run_amplitude_sweep_script) spends ~19k closed-loop
%  runs resolving the (sigma_s, J/N) plane, and every one of them uses the SAME
%  array and the SAME jammer position. Array geometry and jammer angle — the
%  two axes most likely to hide a failure — had never been varied at all. This
%  script varies exactly those, at a fixed nominal amplitude, and scores every
%  case with ONE number so that "is it working" has an answer you can read off
%  a picture instead of arbitrating between nine heatmap panels.
%
%  THE HEADLINE METRIC (plan Section 5, KPI 6)
%      track_score_pct = 100 * mean( (oracle_sinr_db - sinr_db) <= track_tol_db )
%  The fraction of the run spent within track_tol_db of the perfect-knowledge
%  LCMV. It is normalized against the best ACHIEVABLE, which is what makes an
%  easy cell and a hard cell comparable — raw availability saturates at 100%
%  across every easy cell and so cannot rank algorithms there. The oracle
%  scores exactly 100 by construction, so every grid carries its own anchor.
%
%  TWO PROFILES
%      'A'  acceptance grid — the EASY cases, on every array, at every
%           reasonable angle. If this does not pass there is no point reading
%           edge cases, so treat it as the gate on profile B.
%      'B'  stress axes — near-guard, near-endfire, a grating-lobe hunt, low
%           signal, fast drift, and on/off cycle stability. Each axis perturbs
%           profile A along exactly ONE dimension so a failure is attributable.
%
%  READING THE OUTPUT
%      scorecard_<SCN>.png    array x jammer position, score in each cell,
%                             heavy border on anything below the pass mark.
%                             This is the figure to look at first.
%      difficulty_scatter.png score vs steering coherence. Coherence is the
%                             array-INDEPENDENT difficulty coordinate (60 deg
%                             of separation is coherence 0.72 on the 6-element
%                             array and 0.03 on the 20-element one, so raw
%                             angle is not comparable across arrays). Points
%                             below the trend are the blind spots.
%      trace_<SCN>_worst.png  the lowest-scoring case of each scenario,
%                             re-simulated and unrolled in time.
%      acceptance_grid.csv    every metric, seed-mean AND seed-std.
%
%  ON SEEDS
%  Each metric is written twice: <metric> and <metric>_std over the seeds. The
%  amplitude sweep averaged its seeds inline and kept nothing, so the question
%  "are 5 seeds too many" had never been answerable from its output. Read the
%  _std columns and set n_seeds from them rather than from a guess.
%
%  Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P12].
% ══════════════════════════════════════════════════════════════════

% [P12b] Wrapper override hook. A caller (run_mode_c_campaign_script) may set
% a struct GRID_OVERRIDE in the base workspace before running this script to
% vary the suite knobs WITHOUT editing this file — so a campaign that sweeps
% algorithms or adapt config re-uses this exact case-building and scoring code
% rather than a copy of it that can drift. Any field name not in the whitelist
% below is a hard error (CLAUDE.md rule 4: no silent defaults).
if exist('GRID_OVERRIDE', 'var')
    grid_override__ = GRID_OVERRIDE;
else
    grid_override__ = struct();
end
clearvars -except grid_override__; clc;

% ── 0. Paths: reuse the Milestone-1 port + the anti-jam library ─────
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(fileparts(script_dir));           % <repo> = MATLAB/../
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

% ══════════════════════════════════════════════════════════════════
%  SUITE CONFIG — everything you would normally want to edit is here
% ══════════════════════════════════════════════════════════════════

% 'A' = acceptance grid (easy cases), 'B' = stress axes. See the header.
profile = 'A';

% Arrays under test. Polarization and target are declared EXPLICITLY per array
% (CLAUDE.md rule 4: no inherited default).
%
% [P12] ManyDipoles runs 'Theta', NOT 'total' like the other two. The first
% version of this grid held polarization identical across all three for
% comparability, and that was wrong: ManyDipoles is an ideal-dipole export with
% no E_phi, so its Phi component is numerical residue 107 dB down. Pairing it
% with the real Theta component made n_comp = 2 while each source is physically
% rank-1, which broke adapt_music_doa's presence detector (its hardcoded
% n_sig = 2*n_comp then reads presence off a pure-noise eigenvalue) and cost
% the predictive nuller 85% of its detections. select_polarization_stacks now
% rejects that pairing outright.
%
% CAVEAT this creates, and it is not fixable by choosing differently: the
% ManyDipoles row is a single-polarization problem (n_comp = 1) while the other
% two rows are dual-polarization (n_comp = 2). Those are different signal
% models — the SINR definition splits power across components differently — so
% ABSOLUTE scores are not comparable between the ManyDipoles row and the rows
% above it. What stays comparable is the comparison WITHIN a row: one algorithm
% against another, or one jammer position against another.
%
% Grid resolution is NOT uniform either — patchs_with_monopoles and spacing0.6
% export on a 1 deg grid, ManyDipoles on 5 deg. Each array is simulated on its
% own grid, but the angular resolution of a null differs between rows, so the
% preflight prints it and the CSV records it.
array_specs = { ...
    struct('id', 'patchs_with_monopoles', 'dir', 'data/patchs_with_monopoles/', ...
           'polarization', 'total', 'theta_s_deg', 90.0, 'phi_s_deg', 260.0), ...
    struct('id', 'spacing0.6',            'dir', 'data/spacing0.6/', ...
           'polarization', 'total', 'theta_s_deg', 90.0, 'phi_s_deg', 260.0), ...
    struct('id', 'ManyDipoles',           'dir', 'data/ManyDipoles/', ...
           'polarization', 'Theta', 'theta_s_deg', 90.0, 'phi_s_deg', 260.0)};

% Nominal amplitudes for profile A [dB re the sigma_n^2 = 1 noise floor].
% Deliberately comfortable: sigma_s 10 dB clears the 10 dB SINR threshold with
% headroom and sits well clear of the sigma_s <= 4 dB corner P11 identified as
% adaptive loading's weak spot. These are the EASY cases; profile B is where
% they get hard.
nominal_sigma_s_db  = 10.0;
nominal_jn_ratio_db = 20.0;

% Algorithms under test. 'oracle' is always run and is never listed here: it
% defines the reference every score is measured against, so it is not a
% competitor. The first entry is the PRIMARY algorithm — it is the one the
% worst-case trace figures are chosen by.
%
% Cost warning before adding 'predict': it runs a MUSIC eigendecomposition over
% the full far-field grid at every step (adapt.predict.doa_stride = 1), which
% measures 63x lcmv on patchs_with_monopoles and 147x on spacing0.6 (1 deg
% grids, 181x360 points) against 8x on ManyDipoles (5 deg grid). Raise
% adapt.predict.doa_stride to subsample that grid if the full profile is too
% slow — at the cost of DoA resolution.
algorithms = {'lcmv'};

% Restrict the profile to a subset of its scenarios. {} runs all of them.
%
% This exists because algorithm cost is wildly uneven. 'predict' is 63-147x
% lcmv on the 1 deg arrays, and WINDOW runs are 180 s against DRIFT's 30 s, so
% a full two-algorithm profile A costs ~2.1 h of which ~8/9 is spent
% re-confirming that both algorithms score 95-100% on the static scenarios.
% Narrowing to the column that carries the open question answers it in ~14 min.
scenario_filter = {'DRIFT'};

% Oracle-tracking tolerance and the pass mark for the headline score.
track_tol_db = 3.0;
pass_pct     = 90.0;

% Monte Carlo seeds per case. Provisional: with the jammer angle fixed per
% case, the seed varies only the snapshot-noise realization, which should
% average out quickly away from a failure cliff. Read the _std columns in the
% CSV and revise this number from them — that is what they are for.
n_seeds = 3;

% Fraction of the run treated as "steady state" for the _ss metrics.
ss_start_frac = 0.5;

% Route through kpi_evaluate (adds null-pointing error, ~20x slower). The
% inline path in kpi_sweep_metrics mirrors its definitions verbatim.
full_kpi = false;

% Extra output-folder tag, so several campaigns in one session are tellable
% apart on disk. '' = no tag.
out_tag = '';

% Overrides merged into config.adapt for EVERY run of this campaign (e.g.
% forgetting_lambda, or removing loading_factor_db to get fixed loading).
% Fields are merged over the parsed config.yaml block; the special value
% '<remove>' deletes a key instead of setting it, which is how the opt-in
% adaptive loading is turned OFF.
adapt_override = struct();

% ── Apply the wrapper overrides (see the GRID_OVERRIDE hook at the top) ──
if isfield(grid_override__, 'profile'),             profile             = grid_override__.profile;             end
if isfield(grid_override__, 'array_specs'),         array_specs         = grid_override__.array_specs;         end
if isfield(grid_override__, 'algorithms'),          algorithms          = grid_override__.algorithms;          end
if isfield(grid_override__, 'scenario_filter'),     scenario_filter     = grid_override__.scenario_filter;     end
if isfield(grid_override__, 'nominal_sigma_s_db'),  nominal_sigma_s_db  = grid_override__.nominal_sigma_s_db;  end
if isfield(grid_override__, 'nominal_jn_ratio_db'), nominal_jn_ratio_db = grid_override__.nominal_jn_ratio_db; end
if isfield(grid_override__, 'track_tol_db'),        track_tol_db        = grid_override__.track_tol_db;        end
if isfield(grid_override__, 'pass_pct'),            pass_pct            = grid_override__.pass_pct;            end
if isfield(grid_override__, 'n_seeds'),             n_seeds             = grid_override__.n_seeds;             end
if isfield(grid_override__, 'ss_start_frac'),       ss_start_frac       = grid_override__.ss_start_frac;       end
if isfield(grid_override__, 'full_kpi'),            full_kpi            = grid_override__.full_kpi;            end
if isfield(grid_override__, 'out_tag'),             out_tag             = grid_override__.out_tag;             end
if isfield(grid_override__, 'adapt_override'),      adapt_override      = grid_override__.adapt_override;      end
KNOWN_OVERRIDES = {'profile', 'array_specs', 'algorithms', 'scenario_filter', ...
    'nominal_sigma_s_db', 'nominal_jn_ratio_db', 'track_tol_db', 'pass_pct', ...
    'n_seeds', 'ss_start_frac', 'full_kpi', 'out_tag', 'adapt_override'};
unknown__ = setdiff(fieldnames(grid_override__), KNOWN_OVERRIDES);
if ~isempty(unknown__)
    error('run_acceptance_grid:UnknownOverride', ...
        'GRID_OVERRIDE has unrecognised field(s): %s. Known: %s.', ...
        strjoin(unknown__', ', '), strjoin(KNOWN_OVERRIDES, ', '));
end

% ══════════════════════════════════════════════════════════════════

% ── 1. Base config ─────────────────────────────────────────────────
config  = read_config_yaml(fullfile(repo_root, 'config.yaml'));
aj_base = config.antijam;

% ── 2. Build the case list for the selected profile ────────────────
switch upper(profile)
    case 'A'
        cases = build_profile_a(array_specs, nominal_sigma_s_db, nominal_jn_ratio_db);
        profile_name = 'A (acceptance grid — easy cases)';
    case 'B'
        cases = build_profile_b(array_specs, nominal_sigma_s_db, nominal_jn_ratio_db, ...
            aj_base.guard_deg);
        profile_name = 'B (stress axes — edge cases)';
    otherwise
        error('run_acceptance_grid:BadProfile', ...
            'profile must be ''A'' or ''B''; got ''%s''.', profile);
end
if ~isempty(scenario_filter)
    keep = ismember({cases.scn_id}, scenario_filter);
    if ~any(keep)
        error('run_acceptance_grid:EmptyFilter', ...
            ['scenario_filter {%s} matches no scenario in profile %s, whose ' ...
             'scenarios are {%s}.'], strjoin(scenario_filter, ', '), ...
            upper(profile), strjoin(unique({cases.scn_id}, 'stable'), ', '));
    end
    fprintf('Scenario filter: {%s} — %d of %d cases kept.\n', ...
        strjoin(scenario_filter, ', '), sum(keep), numel(cases));
    cases = cases(keep);
end
n_cases = numel(cases);

% ── 3. Output folder ───────────────────────────────────────────────
timestamp  = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
folder_name = sprintf('%s_profile%s', timestamp, upper(profile));
if ~isempty(out_tag)
    folder_name = [folder_name '_' out_tag];
end
output_dir = fullfile(repo_root, 'results', 'acceptance_grid', folder_name);
if ~isfolder(output_dir), mkdir(output_dir); end
fprintf('Profile: %s\n', profile_name);
fprintf('Output:  %s\n', output_dir);

% ── 4. Load each array ONCE and preflight the geometry ─────────────
% Loading patterns per case would dominate the runtime, and the preflight has
% to happen before any simulation: a jammer angle that snaps somewhere
% unintended, or lands inside the guard sector, must stop the run rather than
% quietly produce a plausible-looking number.
fprintf('\nPreflight — loading arrays and resolving requested geometry:\n');
arrays = struct();
for i = 1:numel(array_specs)
    spec = array_specs{i};
    key  = matlab.lang.makeValidName(spec.id);
    patterns = load_element_patterns(fullfile(repo_root, spec.dir));
    cfg_pol  = config;
    cfg_pol.polarization = spec.polarization;
    [s1, s2, pol_label] = select_polarization_stacks(patterns, cfg_pol);

    a = struct();
    a.spec      = spec;
    a.stack1    = s1;
    a.stack2    = s2;
    a.theta_deg = patterns(1).theta_deg;
    a.phi_deg   = patterns(1).phi_deg;
    a.n_el      = size(s1, 1);
    a.pol_label = pol_label;
    a.d_theta   = mean(diff(a.theta_deg));
    a.d_phi     = mean(diff(a.phi_deg));
    a.dir_ref_dbi = NaN;      % filled from this array's first oracle run
    arrays.(key)  = a;

    fprintf('  %-24s %2d el, grid %dx%d (%.0f/%.0f deg step), pol %s, target (%.0f, %.0f)\n', ...
        spec.id, a.n_el, numel(a.theta_deg), numel(a.phi_deg), ...
        a.d_theta, a.d_phi, pol_label, spec.theta_s_deg, spec.phi_s_deg);
end

% Resolve every case's jammer angle now, and record what the ENGINE will
% actually see (nearest_index_2d snaps to the grid) rather than what was asked
% for. Coherence is computed here too, since it depends only on geometry.
fprintf('\nPreflight — jammer geometry per case:\n');
for ic = 1:n_cases
    key = matlab.lang.makeValidName(cases(ic).array_id);
    a   = arrays.(key);
    [th_j, ph_j] = resolve_jammer_angle(cases(ic).sep_deg, cases(ic).cut, ...
        cases(ic).theta_s_deg, cases(ic).phi_s_deg, cases(ic).pos_id);

    % Snap to this array's grid — the same snap sim_engine_init performs — and
    % measure the separation that results. On a 5 deg grid a requested 20 deg
    % can land at 20 exactly; on a coarser one it would not, and reporting the
    % requested value would be a lie.
    [it, ip] = nearest_index_2d(a.theta_deg, a.phi_deg, th_j, ph_j);
    th_snap  = a.theta_deg(it);
    ph_snap  = a.phi_deg(ip);
    sep_true = angular_separation_deg(cases(ic).theta_s_deg, cases(ic).phi_s_deg, ...
        th_snap, ph_snap);

    % Two hard errors, no silent fallbacks. Guard first: the milestone declares
    % a jammer inside the main beam out of scope, so a case that lands there is
    % a mis-specified case, not a hard test.
    if sep_true < aj_base.guard_deg
        error('run_acceptance_grid:AngleInGuard', ...
            ['Case %s / %s snaps to (%.1f, %.1f) deg, only %.2f deg from the ' ...
             'target — inside the %.1f deg guard sector. A jammer in the main ' ...
             'beam is out of scope for this milestone.'], ...
            cases(ic).array_id, cases(ic).pos_id, th_snap, ph_snap, ...
            sep_true, aj_base.guard_deg);
    end
    % Then the snap itself. Half a grid diagonal is the most the nearest-point
    % snap can legitimately move a request; more than that means the requested
    % angle was not representable and the case is measuring something else.
    snap_tol = hypot(a.d_theta, a.d_phi);
    if abs(sep_true - abs(cases(ic).sep_deg)) > snap_tol
        error('run_acceptance_grid:SnapTooFar', ...
            ['Case %s / %s asked for %.1f deg separation but the nearest grid ' ...
             'point on this array gives %.2f deg (grid step %.0f/%.0f deg). ' ...
             'Pick a separation the array grid can represent.'], ...
            cases(ic).array_id, cases(ic).pos_id, abs(cases(ic).sep_deg), ...
            sep_true, a.d_theta, a.d_phi);
    end

    cases(ic).theta_j_deg = th_snap;
    cases(ic).phi_j_deg   = ph_snap;
    cases(ic).sep_true_deg = sep_true;

    % ── Trajectory preflight ──────────────────────────────────────
    % Checking only t = 0 is not enough, and the first run of this suite proved
    % it: a theta-drifting jammer launched 90 deg from the target walks STRAIGHT
    % THROUGH the main beam. sim_scenario's guard_clamp does not complain — it
    % silently pins the offending samples to the guard boundary — so the run
    % completes and reports a plausible-looking availability collapse that is
    % really just "the jammer was sitting on the target for 10% of the run".
    % A jammer inside the main beam is out of scope for this milestone, so a
    % case whose TRAJECTORY enters the guard is a mis-specified case.
    %
    % For a drifting scenario the sign is chosen here rather than in the case
    % builder: it depends on where the jammer starts relative to the target,
    % which is only known after the angle is resolved. Both signs are tried and
    % the one that stays furthest from the target wins.
    [cases(ic).scn, cases(ic).min_sep_deg] = orient_drift( ...
        cases(ic).scn, th_snap, ph_snap, c_aj(aj_base, cases(ic)), config.sim);
    if cases(ic).min_sep_deg < aj_base.guard_deg + 1.0
        error('run_acceptance_grid:TrajectoryEntersGuard', ...
            ['Case %s / %s / %s: the jammer trajectory closes to %.2f deg of the ' ...
             'target (guard is %.1f deg). sim_scenario would clamp those samples ' ...
             'to the guard boundary and the run would silently measure a jammer ' ...
             'sitting in the main beam, which is out of scope. Reduce the drift ' ...
             'rate or duration, or start the jammer further out.'], ...
            cases(ic).array_id, cases(ic).pos_id, cases(ic).scn_id, ...
            cases(ic).min_sep_deg, aj_base.guard_deg);
    end
    cases(ic).coh = kpi_steering_coherence(a.stack1, a.stack2, ...
        a.theta_deg, a.phi_deg, cases(ic).theta_s_deg, cases(ic).phi_s_deg, ...
        th_snap, ph_snap);
    cases(ic).grid_step_deg = a.d_theta;
end

% Print the geometry table once, grouped by array. This table is the answer to
% "why is that row harder than this one" and is worth reading before the
% scorecards.
pos_seen = {};
for ic = 1:n_cases
    tag = [cases(ic).array_id '|' cases(ic).pos_id];
    if any(strcmp(pos_seen, tag)), continue; end
    pos_seen{end + 1} = tag;                                          %#ok<SAGROW>
    fprintf('  %-24s %-12s -> (%6.1f, %6.1f) deg,  sep %6.2f deg,  coherence %.4f\n', ...
        cases(ic).array_id, cases(ic).pos_id, cases(ic).theta_j_deg, ...
        cases(ic).phi_j_deg, cases(ic).sep_true_deg, cases(ic).coh);
end

fprintf('\n%d cases x %d seeds x %d runs (oracle + %s) = %d closed-loop runs\n\n', ...
    n_cases, n_seeds, 1 + numel(algorithms), strjoin(algorithms, ' + '), ...
    n_cases * n_seeds * (1 + numel(algorithms)));

% ── 5. Run ─────────────────────────────────────────────────────────
metric_names = kpi_sweep_metric_names();
csv_path = fullfile(output_dir, 'acceptance_grid.csv');
csv_fid  = fopen(csv_path, 'w');
if csv_fid < 0
    error('run_acceptance_grid:CsvOpen', 'Cannot open %s for writing.', csv_path);
end
fprintf(csv_fid, '%s\n', csv_header(metric_names));

results = repmat(struct('case_index', [], 'algorithm', '', 'mean', [], 'std', []), 0, 1);
n_failed = 0;
t_all = tic;
i_run = 0;
n_runs_total = n_cases * n_seeds * (1 + numel(algorithms));

for ic = 1:n_cases
    c   = cases(ic);
    key = matlab.lang.makeValidName(c.array_id);
    a   = arrays.(key);

    aj = aj_base;
    aj.theta_s_deg = c.theta_s_deg;
    aj.phi_s_deg   = c.phi_s_deg;
    aj.sigma_s_db  = c.sigma_s_db;
    aj.jn_ratio_db = c.jn_ratio_db;

    scn_cfg = c.scn;
    scn_cfg.theta_j_deg = c.theta_j_deg;
    scn_cfg.phi_j_deg   = c.phi_j_deg;
    scn_cfg.jn_ratio_db = c.jn_ratio_db;

    cfg_run = config;
    cfg_run.polarization = a.spec.polarization;
    % [P12b] Campaign-level adapt overrides (forgetting_lambda, loading mode,
    % ...). '<remove>' deletes the key, which is how the opt-in adaptive
    % loading (adapt.loading_factor_db) is switched OFF for a fixed-loading arm.
    ovr_names = fieldnames(adapt_override);
    for iov = 1:numel(ovr_names)
        val = adapt_override.(ovr_names{iov});
        if ischar(val) && strcmp(val, '<remove>')
            if isfield(cfg_run.adapt, ovr_names{iov})
                cfg_run.adapt = rmfield(cfg_run.adapt, ovr_names{iov});
            end
        else
            cfg_run.adapt.(ovr_names{iov}) = val;
        end
    end

    acc = struct('oracle', []);
    for ia = 1:numel(algorithms)
        acc.(algorithms{ia}) = [];
    end
    for iseed = 1:n_seeds
        sim_cfg      = config.sim;
        sim_cfg.seed = config.sim.seed + iseed - 1;
        % One case-seed is wrapped as a unit, matching the amplitude sweep's
        % robustness pattern: a single bad case must not cost the suite, and it
        % must not pass silently either.
        try
            scn = sim_scenario(scn_cfg, aj, sim_cfg);

            o_log = closed_loop_run('oracle', a.stack1, a.stack2, a.theta_deg, ...
                a.phi_deg, scn, aj, sim_cfg, cfg_run, []);
            o_log.oracle_sinr_db = o_log.sinr_db;
            i_run = i_run + 1;

            % Quiescent directivity is a property of the ARRAY, so it is
            % computed once per array — not once per campaign, which is what
            % the amplitude sweep did and what silently breaks the moment a
            % second array is in play (5.13 / 4.61 / 17.16 dBi across these
            % three, i.e. a 12 dB error on the wrong reference).
            if ~isfinite(arrays.(key).dir_ref_dbi)
                arrays.(key).dir_ref_dbi = kpi_quiescent_directivity(o_log, ...
                    a.n_el, a.stack1, a.stack2, a.theta_deg, a.phi_deg);
                a.dir_ref_dbi = arrays.(key).dir_ref_dbi;
                fprintf('  [%s] quiescent-beam directivity toward target: %.2f dBi\n', ...
                    c.array_id, a.dir_ref_dbi);
            end

            acc.oracle = accumulate_metrics(acc.oracle, metric_names, ...
                kpi_sweep_metrics(o_log, scn, aj, a.stack1, a.stack2, ...
                    a.theta_deg, a.phi_deg, ss_start_frac, full_kpi, ...
                    a.dir_ref_dbi, track_tol_db));

            % Every algorithm sees the SAME scenario and the same seed, and is
            % scored against the same oracle run — so a difference between two
            % rows of the scorecard is a difference between the algorithms and
            % not between the channels they happened to be handed.
            for ia = 1:numel(algorithms)
                g_log = closed_loop_run(algorithms{ia}, a.stack1, a.stack2, ...
                    a.theta_deg, a.phi_deg, scn, aj, sim_cfg, cfg_run, []);
                g_log.oracle_sinr_db = o_log.sinr_db;
                i_run = i_run + 1;
                acc.(algorithms{ia}) = accumulate_metrics( ...
                    acc.(algorithms{ia}), metric_names, ...
                    kpi_sweep_metrics(g_log, scn, aj, a.stack1, a.stack2, ...
                        a.theta_deg, a.phi_deg, ss_start_frac, full_kpi, ...
                        a.dir_ref_dbi, track_tol_db));
            end
        catch err
            n_failed = n_failed + 1;
            warning('run_acceptance_grid:CaseFailed', ...
                ['Case %s / %s / %s, seed %d FAILED: %s (%s). That seed is ' ...
                 'dropped; the case keeps its other seeds.'], ...
                c.array_id, c.pos_id, c.scn_id, sim_cfg.seed, ...
                err.message, err.identifier);
        end
    end

    score_txt = '';
    for alg = [{'oracle'}, algorithms]
        [mu, sd] = seed_stats(acc.(alg{1}), metric_names);
        results(end + 1) = struct('case_index', ic, 'algorithm', alg{1}, ...
            'mean', mu, 'std', sd);                                   %#ok<SAGROW>
        fprintf(csv_fid, '%s\n', csv_row(c, alg{1}, mu, sd, metric_names, a.n_el));
        if ~strcmp(alg{1}, 'oracle')
            score_txt = [score_txt sprintf(' %s %5.1f%%', alg{1}, ...
                mu.track_score_pct)];                                 %#ok<AGROW>
        end
    end

    elapsed = toc(t_all);
    fprintf('  [%2d/%2d] %-22s %-12s %-14s %s  (%d/%d runs, %.0f s, ~%.0f s left)\n', ...
        ic, n_cases, c.array_id, c.pos_id, c.scn_id, score_txt, ...
        i_run, n_runs_total, elapsed, elapsed * (n_runs_total - i_run) / max(i_run, 1));
end
fclose(csv_fid);
fprintf('\nSuite complete in %.0f s (%d failed case-seeds).\n', toc(t_all), n_failed);

grid_out = struct();
grid_out.profile      = upper(profile);
grid_out.cases        = cases;
grid_out.results      = results;
grid_out.metric_names = metric_names;
grid_out.algorithms   = algorithms;
grid_out.n_seeds      = n_seeds;
grid_out.track_tol_db = track_tol_db;
grid_out.pass_pct     = pass_pct;
grid_out.ss_start_frac = ss_start_frac;
grid_out.adapt_override = adapt_override;      % [P12b] which arm this was
grid_out.out_tag        = out_tag;
grid_out.nominal_sigma_s_db  = nominal_sigma_s_db;
grid_out.nominal_jn_ratio_db = nominal_jn_ratio_db;
save(fullfile(output_dir, 'acceptance_grid.mat'), 'grid_out');

% ── 6. Figures ─────────────────────────────────────────────────────
% Data is on disk before any of this runs, so a graphics failure costs a
% picture and not the campaign.
fprintf('\nRendering figures...\n');
scenario_ids = unique({cases.scn_id}, 'stable');
array_ids    = unique({cases.array_id}, 'stable');
pos_ids      = unique({cases.pos_id}, 'stable');

for is = 1:numel(scenario_ids)
    % One panel per algorithm on a shared colour scale, so two algorithms are
    % compared by looking across the figure rather than by flipping between
    % two files.
    panels = struct('map', {}, 'title', {}, 'note', {});
    for ia = 1:numel(algorithms)
        map  = NaN(numel(array_ids), numel(pos_ids));
        note = repmat({''}, numel(array_ids), numel(pos_ids));
        for ir = 1:numel(results)
            if ~strcmp(results(ir).algorithm, algorithms{ia}), continue; end
            c = cases(results(ir).case_index);
            if ~strcmp(c.scn_id, scenario_ids{is}), continue; end
            iy = find(strcmp(array_ids, c.array_id), 1);
            ix = find(strcmp(pos_ids,   c.pos_id),   1);
            map(iy, ix)  = results(ir).mean.track_score_pct;
            note{iy, ix} = sprintf('c=%.2f', c.coh);
        end
        % Built by field assignment rather than struct(...): struct() consumes
        % a cell argument to define a struct ARRAY, so 'note', {note} would
        % produce one panel per annotation cell and 'note', {{note}} would nest
        % the grid one level too deep. Neither is what is wanted here.
        p = struct();
        p.map   = map;
        p.title = sprintf('%s — %s', scenario_ids{is}, algorithms{ia});
        p.note  = note;
        panels(ia) = p;
    end
    safe_plot(@() plot_scorecard(array_ids, pos_ids, panels, sprintf( ...
        ['Oracle-tracking score — profile %s, scenario %s ' ...
         '(%d seeds, sigma_s %.0f dB, J/N %.0f dB; c = steering coherence)'], ...
        upper(profile), scenario_ids{is}, n_seeds, nominal_sigma_s_db, ...
        nominal_jn_ratio_db), ...
        fullfile(output_dir, sprintf('scorecard_%s.png', scenario_ids{is})), pass_pct));
end

for ia = 1:numel(algorithms)
    points = struct('coh', {}, 'score_pct', {}, 'sep_deg', {}, 'array_id', {}, ...
                    'scenario_id', {}, 'label', {});
    for ir = 1:numel(results)
        if ~strcmp(results(ir).algorithm, algorithms{ia}), continue; end
        c = cases(results(ir).case_index);
        points(end + 1) = struct('coh', c.coh, ...
            'score_pct', results(ir).mean.track_score_pct, ...
            'sep_deg', c.sep_true_deg, 'array_id', c.array_id, ...
            'scenario_id', c.scn_id, ...
            'label', sprintf('%s/%s/%s', c.array_id, c.pos_id, c.scn_id)); %#ok<SAGROW>
    end
    safe_plot(@() plot_difficulty_scatter(points, sprintf( ...
        'Blind-spot map — profile %s (%s, %d seeds)', upper(profile), ...
        algorithms{ia}, n_seeds), ...
        fullfile(output_dir, sprintf('difficulty_scatter_%s.png', ...
            algorithms{ia})), pass_pct));
end

% Worst case per scenario, re-simulated at the base seed and unrolled in time.
% The scorecard says WHICH case is bad; only a trace says why.
for is = 1:numel(scenario_ids)
    % Picked on the primary algorithm (the first entry of `algorithms`); the
    % trace itself draws every algorithm, so the comparison is still visible.
    sel = find(arrayfun(@(r) strcmp(r.algorithm, algorithms{1}) && ...
        strcmp(cases(r.case_index).scn_id, scenario_ids{is}), results));
    if isempty(sel), continue; end
    scores = arrayfun(@(k) results(k).mean.track_score_pct, sel);
    scores(~isfinite(scores)) = -Inf;      % a failed case is the worst case
    [~, k] = min(scores);
    safe_plot(@() render_worst_trace(cases(results(sel(k)).case_index), ...
        arrays, config, aj_base, output_dir, scenario_ids{is}, algorithms));
end

fprintf('\nDone. All artifacts in:\n  %s\n', output_dir);
fprintf(['\nRead in this order:\n' ...
         '  1. scorecard_<SCENARIO>.png   — which (array, angle) cells pass\n' ...
         '  2. difficulty_scatter.png     — which failures are real blind spots\n' ...
         '  3. trace_<SCENARIO>_worst.png — what the worst case actually did\n' ...
         '  4. acceptance_grid.csv        — the _std columns set n_seeds\n']);


% ────────────────────────── CASE BUILDERS ─────────────────────────

function cases = build_profile_a(array_specs, sigma_s_db, jn_ratio_db)
% Profile A: the easy cases. One nominal amplitude, three benign scenarios,
% and jammer positions specified as SEPARATION from the target rather than as
% absolute angles — so the same five positions mean the same thing on every
% array, which is what makes the scorecard rows comparable.
%
% sep45_ph is the same 45 deg separation as sep45_th, displaced into the
% orthogonal cut. It is there to catch axis-dependent behaviour: an array whose
% theta and phi cuts differ will score them differently, and a single-cut
% position set would never notice.
%
% sep135 has to be a phi-cut position. With the target at theta_s = 90 the
% theta cut spans only +-90 deg of separation before running off the [0, 180]
% elevation range, so a 135 deg separation is unreachable there — the preflight
% raises run_acceptance_grid:UnreachableSeparation rather than silently
% clamping it to 90. The phi cut reaches 180 deg at theta_s = 90.
positions = { ...
    struct('id', 'sep20_th',  'sep_deg',  20.0, 'cut', 'theta'), ...
    struct('id', 'sep45_th',  'sep_deg',  45.0, 'cut', 'theta'), ...
    struct('id', 'sep90_th',  'sep_deg',  90.0, 'cut', 'theta'), ...
    struct('id', 'sep135_ph', 'sep_deg', 135.0, 'cut', 'phi'), ...
    struct('id', 'sep45_ph',  'sep_deg',  45.0, 'cut', 'phi')};

% STATIC is the baseline. DRIFT is the milestone's primary tracking scenario.
% WINDOW is the clean lifecycle measurement: one turn_on and one turn_off with
% quiet either side, so recovery AND release are both visible without the
% averaging that a multi-cycle scenario forces. Multi-cycle ONOFF lives in
% profile B, where its per-cycle behaviour is the actual question.
scns = { ...
    struct('id', 'STATIC', 'motion', 'static', 'power', 'constant', ...
           'duration_s', 60.0), ...
    ... % DRIFT runs 30 s, not the 60 s the amplitude sweep used: at 2 deg/s
    ... % that is 60 deg of jammer travel, which is enough to test tracking and
    ... % short enough that no position's trajectory reaches the guard sector.
    ... % A 60 s drift launched from sep90_th (theta_j = 180, the fold point)
    ... % walks the jammer straight through the main beam whichever way it is
    ... % pointed, and the trajectory preflight rejects it.
    struct('id', 'DRIFT',  'motion', 'drift',  'power', 'constant', ...
           'theta_drift_deg_per_s', 2.0, 'phi_drift_deg_per_s', 0.0, ...
           'duration_s', 30.0), ...
    struct('id', 'WINDOW', 'motion', 'static', 'power', 'window', ...
           'on_time_s', 60.0, 'off_time_s', 120.0, 'duration_s', 180.0)};

cases = cross_cases(array_specs, positions, scns, sigma_s_db, jn_ratio_db);
end


function cases = build_profile_b(array_specs, sigma_s_db, jn_ratio_db, guard_deg)
% Profile B: stress axes. Each block perturbs profile A along exactly ONE
% dimension, so a failure points at a cause instead of at a combination.
cases = repmat(empty_case(), 0, 1);

scn_static = struct('id', 'STATIC', 'motion', 'static', 'power', 'constant', ...
                    'duration_s', 60.0);
scn_window = struct('id', 'WINDOW', 'motion', 'static', 'power', 'window', ...
                    'on_time_s', 60.0, 'off_time_s', 120.0, 'duration_s', 180.0);

% (1) Near-guard and near-endfire. The closest the jammer is allowed to get,
%     and the two directions where element patterns roll off. With the target
%     at theta_s = 90, a theta-cut separation of -85 / +85 deg puts the jammer
%     at theta 5 and theta 175 respectively.
edge_pos = { ...
    struct('id', 'near_guard', 'sep_deg', guard_deg + 5.0, 'cut', 'theta'), ...
    struct('id', 'endfire_lo', 'sep_deg', -85.0,           'cut', 'theta'), ...
    struct('id', 'endfire_hi', 'sep_deg',  85.0,           'cut', 'theta')};
cases = [cases; cross_cases(array_specs, edge_pos, {scn_static, scn_window}, ...
    sigma_s_db, jn_ratio_db)];

% (2) Grating-lobe hunt, spacing0.6 only. A 4x4 array at 0.6 lambda should fold
%     some angles back toward the target; that shows up as a COHERENCE spike at
%     a large separation, which the difficulty scatter's lower panel plots
%     directly. Fine angular steps here are the whole point.
sp = array_specs(strcmp(cellfun(@(s) s.id, array_specs, 'UniformOutput', false), ...
    'spacing0.6'));
if ~isempty(sp)
    sweep_pos = cell(1, 0);
    for s = 10:5:70
        sweep_pos{end + 1} = struct('id', sprintf('gl%03d', s), ...
            'sep_deg', double(s), 'cut', 'theta');                    %#ok<AGROW>
    end
    cases = [cases; cross_cases(sp, sweep_pos, {scn_static}, sigma_s_db, jn_ratio_db)];
end

% (3) Low signal. sigma_s = 0 dB is the corner P11 measured adaptive loading
%     losing badly in (WINDOW recovery 486 steps vs 12.8 for fixed loading).
low_pos = { ...
    struct('id', 'sep45_th', 'sep_deg', 45.0, 'cut', 'theta')};
low = cross_cases(array_specs, low_pos, {scn_static, scn_window}, 0.0, jn_ratio_db);
for i = 1:numel(low), low(i).scn_id = [low(i).scn_id '_LOWSIG']; end
cases = [cases; low];

% (4) Fast drift — faster than the covariance forgetting time constant, where a
%     tracker that keeps up at 2 deg/s has no reason to keep up at 20.
%
%     This block runs on a PHI-cut position, unlike everything else here. A
%     theta-cut jammer drifting at 10-20 deg/s covers 600-1200 deg in a 60 s
%     run, folding through the main beam over and over; the trajectory
%     preflight rejects it, and rightly so — the measurement would be of a
%     jammer repeatedly transiting the target, not of fast tracking. Started at
%     theta_j = theta_s with the azimuth held 45 deg away, the same theta drift
%     sweeps the jammer over the pole and back with the separation bounded in
%     [45, 90] deg, so the motion is fast and the geometry never degenerates.
drift_pos = { ...
    struct('id', 'sep45_ph', 'sep_deg', 45.0, 'cut', 'phi')};
for rate = [10.0, 20.0]
    scn_fast = struct('id', sprintf('DRIFT%02d', rate), 'motion', 'drift', ...
        'power', 'constant', 'theta_drift_deg_per_s', rate, ...
        'phi_drift_deg_per_s', 0.0, 'duration_s', 60.0);
    cases = [cases; cross_cases(array_specs, drift_pos, {scn_fast}, ...
        sigma_s_db, jn_ratio_db)];                                    %#ok<AGROW>
end

% (5) Cycle stability. Three on/off cycles, which is the ONLY reason a
%     multi-cycle scenario is kept anywhere in this suite: not to average
%     recovery over more events (seeds do that better) but to check the tracker
%     does not accumulate state across cycles. Read it per cycle in the trace
%     figure — an averaged recovery time is exactly what would hide a drift.
scn_cycles = struct('id', 'ONOFF3', 'motion', 'static', 'power', 'onoff', ...
    'duty_cycle', 0.5, 'toggle_period_s', 20.0, 'duration_s', 120.0);
cases = [cases; cross_cases(array_specs, low_pos, {scn_cycles}, ...
    sigma_s_db, jn_ratio_db)];
end


function cases = cross_cases(array_specs, positions, scns, sigma_s_db, jn_ratio_db)
% Full cross product of arrays x positions x scenarios, as a struct array.
cases = repmat(empty_case(), 0, 1);
for ia = 1:numel(array_specs)
    for ip = 1:numel(positions)
        for is = 1:numel(scns)
            c = empty_case();
            c.array_id    = array_specs{ia}.id;
            c.theta_s_deg = array_specs{ia}.theta_s_deg;
            c.phi_s_deg   = array_specs{ia}.phi_s_deg;
            c.pos_id      = positions{ip}.id;
            c.sep_deg     = positions{ip}.sep_deg;
            c.cut         = positions{ip}.cut;
            c.scn         = scns{is};
            c.scn_id      = scns{is}.id;
            c.sigma_s_db  = sigma_s_db;
            c.jn_ratio_db = jn_ratio_db;
            cases(end + 1, 1) = c;                                    %#ok<AGROW>
        end
    end
end
end


function c = empty_case()
% One blank case. Declared in one place so every builder produces struct arrays
% with identical field sets — MATLAB refuses to concatenate ones that differ.
c = struct('array_id', '', 'theta_s_deg', NaN, 'phi_s_deg', NaN, ...
           'pos_id', '', 'sep_deg', NaN, 'cut', '', ...
           'scn', struct(), 'scn_id', '', ...
           'sigma_s_db', NaN, 'jn_ratio_db', NaN, ...
           'theta_j_deg', NaN, 'phi_j_deg', NaN, 'sep_true_deg', NaN, ...
           'coh', NaN, 'grid_step_deg', NaN, 'min_sep_deg', NaN);
end


% ────────────────────────── GEOMETRY ──────────────────────────────

function [theta_j_deg, phi_j_deg] = resolve_jammer_angle(sep_deg, cut, ...
    theta_s_deg, phi_s_deg, pos_id)
% Turn a requested angular SEPARATION into an absolute (theta, phi) jammer
% direction. Positions are specified as separations so that one position list
% means the same thing on every array; this is where that is made concrete.
%
% 'theta' cut: signed, so a caller can address both sides of the target
%              (sep = -85 with theta_s = 90 gives theta_j = 5, near one
%              endfire; sep = +85 gives 175, near the other).
% 'phi'   cut: theta held at theta_s and the azimuth offset solved from the
%              spherical law of cosines,
%                  cos(sep) = cos^2(theta_s) + sin^2(theta_s) cos(dphi).
%              Not every separation is reachable in this cut — the maximum is
%              2*theta_s — so an unreachable request is an error, not a clamp.
switch lower(cut)
    case 'theta'
        theta_j_deg = theta_s_deg + sep_deg;
        if theta_j_deg > 180 || theta_j_deg < 0
            theta_j_deg = theta_s_deg - sep_deg;
        end
        if theta_j_deg > 180 || theta_j_deg < 0
            error('run_acceptance_grid:UnreachableSeparation', ...
                ['Position ''%s'': a %.1f deg theta-cut separation does not fit ' ...
                 'either side of theta_s = %.1f deg within [0, 180].'], ...
                pos_id, sep_deg, theta_s_deg);
        end
        phi_j_deg = phi_s_deg;
    case 'phi'
        ts = deg2rad(theta_s_deg);
        if abs(sin(ts)) < 1e-9
            error('run_acceptance_grid:PhiCutAtPole', ...
                ['Position ''%s'': a phi-cut separation is undefined at ' ...
                 'theta_s = %.1f deg (every azimuth is the same direction).'], ...
                pos_id, theta_s_deg);
        end
        cos_dphi = (cos(deg2rad(abs(sep_deg))) - cos(ts)^2) / sin(ts)^2;
        if cos_dphi < -1 || cos_dphi > 1
            error('run_acceptance_grid:UnreachableSeparation', ...
                ['Position ''%s'': a %.1f deg separation is not reachable in the ' ...
                 'phi cut at theta_s = %.1f deg (maximum is %.1f deg).'], ...
                pos_id, abs(sep_deg), theta_s_deg, 2 * theta_s_deg);
        end
        theta_j_deg = theta_s_deg;
        phi_j_deg   = mod(phi_s_deg + sign_or_one(sep_deg) * ...
            rad2deg(acos(cos_dphi)), 360);
    otherwise
        error('run_acceptance_grid:BadCut', ...
            'Position ''%s'': cut must be ''theta'' or ''phi''; got ''%s''.', ...
            pos_id, cut);
end
end


function aj = c_aj(aj_base, c)
% The antijam section a case will actually run with — needed by the trajectory
% preflight, which has to build the scenario exactly as the run loop will.
aj = aj_base;
aj.theta_s_deg = c.theta_s_deg;
aj.phi_s_deg   = c.phi_s_deg;
aj.sigma_s_db  = c.sigma_s_db;
aj.jn_ratio_db = c.jn_ratio_db;
end


function [scn_cfg, min_sep_deg] = orient_drift(scn_cfg, theta_j_deg, phi_j_deg, ...
                                               aj, sim_cfg)
% Point a drifting jammer AWAY from the target, and report how close its
% trajectory comes.
%
% Both drift signs are simulated and the one whose minimum separation is
% largest is kept. Picking the sign analytically is tempting but wrong: theta
% folds at 0 and 180 (sim_scenario's reflect_into), so a jammer launched at
% theta = 180 moves back toward the target whichever way it is nudged, and only
% actually evaluating the trajectory settles which nudge is less bad.
%
% Scenarios without a drift are passed through unchanged; they still get their
% minimum separation measured, because a static jammer's separation is a fact
% worth recording next to the score.
scn_cfg.theta_j_deg = theta_j_deg;
scn_cfg.phi_j_deg   = phi_j_deg;
if ~isfield(scn_cfg, 'jn_ratio_db')
    scn_cfg.jn_ratio_db = aj.jn_ratio_db;
end

has_drift = isfield(scn_cfg, 'theta_drift_deg_per_s') && ...
            scn_cfg.theta_drift_deg_per_s ~= 0;
if ~has_drift
    min_sep_deg = min(trajectory_separation(scn_cfg, aj, sim_cfg));
    return
end

rate  = abs(scn_cfg.theta_drift_deg_per_s);
best  = -Inf;
for s = [1, -1]
    trial = scn_cfg;
    trial.theta_drift_deg_per_s = s * rate;
    m = min(trajectory_separation(trial, aj, sim_cfg));
    if m > best
        best    = m;
        scn_cfg = trial;
    end
end
min_sep_deg = best;
end


function sep = trajectory_separation(scn_cfg, aj, sim_cfg)
% Angular separation target-to-jammer at every step of a scenario.
scn = sim_scenario(scn_cfg, aj, sim_cfg);
sep = angular_separation_deg(scn.theta_j_deg, scn.phi_j_deg, ...
    aj.theta_s_deg, aj.phi_s_deg);
end


function s = sign_or_one(v)
% sign(), except sign(0) = 1 — a zero separation still needs a direction.
s = sign(v);
if s == 0, s = 1; end
end


% ────────────────────────── METRIC BOOKKEEPING ────────────────────

function acc = accumulate_metrics(acc, names, m)
% Append one seed's metric struct to a struct of per-seed vectors.
if isempty(acc)
    acc = struct();
    for i = 1:numel(names), acc.(names{i}) = []; end
end
for i = 1:numel(names)
    acc.(names{i})(end + 1) = m.(names{i});
end
end


function [mu, sd] = seed_stats(acc, names)
% Collapse per-seed vectors to a mean AND a standard deviation.
%
% The std is the point. The amplitude sweep averaged its seeds inline and kept
% nothing else, so "is this difference between two loading modes real, or is it
% seed noise" and "are 5 seeds more than we need" were both unanswerable from
% its output. Writing sd beside mu makes both answerable from the CSV alone.
mu = struct(); sd = struct();
for i = 1:numel(names)
    if isempty(acc)
        mu.(names{i}) = NaN;      % every seed of this case failed
        sd.(names{i}) = NaN;
        continue
    end
    v = acc.(names{i});
    v = v(isfinite(v));
    if isempty(v)
        mu.(names{i}) = NaN;
        sd.(names{i}) = NaN;
    else
        mu.(names{i}) = mean(v);
        % std of one sample is NaN in MATLAB; report 0 so a single-seed run
        % does not look like a failed measurement.
        if numel(v) < 2, sd.(names{i}) = 0; else, sd.(names{i}) = std(v); end
    end
end
end


function h = csv_header(names)
cols = {'array_id', 'n_el', 'pos_id', 'sep_deg', 'min_sep_deg', 'coh', ...
        'theta_j_deg', 'phi_j_deg', 'grid_step_deg', 'scenario', 'algorithm', ...
        'sigma_s_db', 'jn_ratio_db', 'js_db'};
for i = 1:numel(names)
    cols{end + 1} = names{i};                                         %#ok<AGROW>
    cols{end + 1} = [names{i} '_std'];                                %#ok<AGROW>
end
h = strjoin(cols, ',');
end


function row = csv_row(c, alg, mu, sd, names, n_el)
vals = cell(1, 2 * numel(names));
for i = 1:numel(names)
    vals{2 * i - 1} = sprintf('%.4g', mu.(names{i}));
    vals{2 * i}     = sprintf('%.4g', sd.(names{i}));
end
row = sprintf('%s,%d,%s,%.2f,%.2f,%.4f,%.1f,%.1f,%g,%s,%s,%g,%g,%g,%s', ...
    c.array_id, n_el, c.pos_id, c.sep_true_deg, c.min_sep_deg, c.coh, c.theta_j_deg, ...
    c.phi_j_deg, c.grid_step_deg, c.scn_id, alg, c.sigma_s_db, ...
    c.jn_ratio_db, c.jn_ratio_db - c.sigma_s_db, strjoin(vals, ','));
end


% ────────────────────────── FIGURES ───────────────────────────────

function render_worst_trace(c, arrays, config, aj_base, ...
                            output_dir, scenario_id, algorithms)
% Re-simulate one case at the base seed and draw its time histories. The
% scorecard says which case is bad; only a trace says why — one long outage and
% a comb of post-toggle dips produce the same scalar.
key = matlab.lang.makeValidName(c.array_id);
a   = arrays.(key);

aj = aj_base;
aj.theta_s_deg = c.theta_s_deg;
aj.phi_s_deg   = c.phi_s_deg;
aj.sigma_s_db  = c.sigma_s_db;
aj.jn_ratio_db = c.jn_ratio_db;

scn_cfg = c.scn;
scn_cfg.theta_j_deg = c.theta_j_deg;
scn_cfg.phi_j_deg   = c.phi_j_deg;
scn_cfg.jn_ratio_db = c.jn_ratio_db;

cfg_run = config;
cfg_run.polarization = a.spec.polarization;
scn = sim_scenario(scn_cfg, aj, config.sim);

o_log = closed_loop_run('oracle', a.stack1, a.stack2, a.theta_deg, a.phi_deg, ...
    scn, aj, config.sim, cfg_run, []);

traces = struct('label', {}, 'sinr_db', {}, 'dir_s_dbi', {}, 'is_oracle', {});
traces(1) = struct('label', 'oracle', 'sinr_db', o_log.sinr_db, ...
    'dir_s_dbi', compute_directivity_trace(o_log, a.stack1, a.stack2, ...
        a.theta_deg, a.phi_deg), ...
    'is_oracle', true);
for ia = 1:numel(algorithms)
    g_log = closed_loop_run(algorithms{ia}, a.stack1, a.stack2, a.theta_deg, ...
        a.phi_deg, scn, aj, config.sim, cfg_run, []);
    traces(end + 1) = struct('label', algorithms{ia}, ...
        'sinr_db', g_log.sinr_db, ...
        'dir_s_dbi', compute_directivity_trace(g_log, a.stack1, a.stack2, ...
            a.theta_deg, a.phi_deg), ...
        'is_oracle', false);                                          %#ok<AGROW>
end

plot_cell_traces(scn, aj, traces, a.dir_ref_dbi, sprintf( ...
    'Worst case for %s — %s / %s (sep %.1f deg, coherence %.3f)', ...
    scenario_id, c.array_id, c.pos_id, c.sep_true_deg, c.coh), ...
    fullfile(output_dir, sprintf('trace_%s_worst.png', scenario_id)));
end


function safe_plot(fn)
% Run one figure call, converting a failure into a warning. The data is already
% on disk by the time any of these run, and a graphics hiccup must not take the
% other figures down with it.
try
    fn();
catch err
    warning('run_acceptance_grid:PlotFailed', ...
        'A figure failed to render (%s): %s', err.identifier, err.message);
end
end
