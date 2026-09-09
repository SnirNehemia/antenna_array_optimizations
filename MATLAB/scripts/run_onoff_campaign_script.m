% ==================================================================
%  run_onoff_campaign_script.m -- [O] stationary on/off jammer, every array
%
%  WHAT THIS ANSWERS
%  How close does each Mode C algorithm get to what a GIVEN array can actually
%  do, against a stationary jammer that switches on and off -- across every
%  array in data/, several target directions, several jammer separations and
%  several toggle periods.
%
%  WHY IT IS BUILT DIFFERENTLY FROM THE ACCEPTANCE GRID
%  Two axes that every previous campaign held fixed are swept here, and both
%  turn out to matter more than the algorithm choice:
%
%    * TARGET DIRECTION. Every campaign to date fixed the target at
%      (theta 90, phi 260). Measured across the seven arrays, the quiescent
%      directivity toward a target swings ~18 dB and the beamwidth spans
%      24-120 deg, so target direction is a difficulty axis in its own right.
%      On ManyDipoles, (90, 260) is additionally the mirror-symmetry plane --
%      that array's single most favourable target.
%
%    * SEPARATION IN UNITS OF THE ARRAY'S OWN BEAM. A jammer 20 deg from the
%      target is trivial for a 20-element array and inside the main lobe of a
%      6-element one. Separations are therefore specified as MULTIPLES of the
%      derived guard (half the wider HPBW toward that target), so the same
%      three numbers mean the same physical difficulty everywhere.
%
%  FEASIBILITY IS MEASURED, NOT ASSUMED
%  `kpi_array_profile` decides per (array, target) whether the cell is a real
%  test at all. Two ways it is not: the array has no gain toward that direction
%  (Monopoles reads -2.45 dBi at (90, 30)), or its main beam covers so much of
%  the sphere that no jammer position lies outside it (every Dipole target).
%  Those cells are recorded with a reason rather than scored, which is what
%  lets a 1-element array be reported honestly instead of as a failure.
%
%  ARMS
%    base   -- the shipped predict path.
%    onoff  -- the [O] repair: presence from a short-memory covariance, an
%              analysis window sized for the longest period worth detecting,
%              and a lead scaled to the detected period and capped at a small
%              number of covariance horizons.
%
%  Part of: Antenna Array Pattern Optimization Tool -- anti-jam milestone [O].
% ==================================================================

function run_onoff_campaign_script(arm_filter)
% ARM_FILTER : optional cellstr of arm tags ({} or absent = all).

if nargin < 1, arm_filter = {}; end
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(fileparts(script_dir));
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

% ══════════════════════════════════════════════════════════════════
%  SUITE CONFIG
% ══════════════════════════════════════════════════════════════════
ARRAYS = {'patchs_with_monopoles', 'spacing0.6', 'spacing0.6_disturbed3', ...
          'Monopoles', 'ManyDipoles', 'patch_back2back', 'Dipole'};

% Target directions, common across arrays so rows are comparable.
% [O2] Trimmed from five to three. (90, 30) and (20, 260) contributed most of
% the infeasible cells and no axis the other three do not cover, and the run
% length below doubled -- the budget was better spent on cycles than targets.
TARGETS = [90 260; 60 260; 45 150];

% Jammer separation as a MULTIPLE of that cell's derived guard.
SEP_MULT = [1.25, 2.5];

% Toggle periods [s]. Chosen to straddle the two failure modes measured in O0:
% below the covariance horizon's ability to resolve OFF phases, and above the
% window length that made a period detectable at all.
PERIODS  = [4.0, 10.0, 25.0];
DUTY     = 0.5;
% [O2] EIGHT cycles, not four. Period detection needs min_periods = 3 cycles,
% so a 4-cycle run spends ~75% of itself unlearned -- which understated the
% anticipatory benefit and, more importantly, stopped the off-window-adaptive
% release from ever engaging (it needs a learned period to size its ramp).
% A fielded jammer toggles for far longer than four cycles; eight is closer to
% the regime the system actually occupies and is what the budget allows.
CYCLES   = 8;

SIGMA_S_DB  = 10.0;
JN_RATIO_DB = 25.0;
N_SEEDS     = 3;
TRACK_TOL_DB = 3.0;   % "within this of the oracle" defines the score

% A cell is a real test only if the array has usable gain toward the target AND
% some direction lies outside its main beam.
MIN_DIR_DBI   = 0.0;
MAX_GUARD_DEG = 90.0;

% The repair as measured in the first campaign: presence from a fast covariance,
% a window sized for the longest period, a period-scaled lead -- and a BINARY
% release of the null the moment the jammer is declared absent.
ONOFF_BLOCK = struct('fast_lambda', 0.5, 'max_period_s', 60.0, ...
                     'lead_frac', 0.15, 'lead_cap_horizons', 2.0);

% [O2] The same, with the release GRADED instead of binary: the null is relaxed
% continuously through diagonal loading as measured absence accumulates, and
% once the period is learned the ramp is sized to complete within a fraction of
% the predicted OFF window. That last part is what lets one setting serve both
% regimes -- the ramp is short relative to a long gap and never completes inside
% a short one, which no fixed threshold can do.
GRADED_BLOCK = ONOFF_BLOCK;
GRADED_BLOCK.release_min_horizons  = 2.0;
GRADED_BLOCK.release_ramp_horizons = 2.0;
GRADED_BLOCK.release_off_frac      = 0.30;

ARMS = { struct('tag', 'base',   'onoff', []), ...
         struct('tag', 'onoff',  'onoff', ONOFF_BLOCK), ...
         struct('tag', 'graded', 'onoff', GRADED_BLOCK) };

tags = cellfun(@(a) a.tag, ARMS, 'UniformOutput', false);
if ~isempty(arm_filter)
    keep = ismember(tags, arm_filter);
    if ~any(keep)
        error('run_onoff_campaign:EmptyFilter', ...
            'arm_filter {%s} matches no arm; known arms are {%s}.', ...
            strjoin(arm_filter, ', '), strjoin(tags, ', '));
    end
    ARMS = ARMS(keep);
end

config = read_config_yaml(fullfile(repo_root, 'config.yaml'));

% ══════════════════════════════════════════════════════════════════
%  1. Load every array once and profile it against every target
% ══════════════════════════════════════════════════════════════════
fprintf('Preflight -- loading arrays and profiling each (array, target):\n\n');
fprintf('  %-22s %-10s %5s %8s %8s %8s %7s %6s  %s\n', ...
    'array', 'target', 'n_el', 'dir dBi', 'HPBW th', 'GUARD', 'mirror', 'MUSIC', 'verdict');

arrays = struct();
cases  = struct('array_id', {}, 'key', {}, 'target', {}, 'prof', {}, ...
                'sep_mult', {}, 'sep_deg', {}, 'theta_j_deg', {}, 'phi_j_deg', {}, ...
                'period_s', {}, 'coh', {});
skipped = {};

for ia = 1:numel(ARRAYS)
    A   = ARRAYS{ia};
    key = matlab.lang.makeValidName(A);
    patterns = load_element_patterns(fullfile(repo_root, 'data', A, filesep));
    pol = pick_polarization(patterns);
    cfg_pol = config; cfg_pol.polarization = pol;
    [s1, s2] = select_polarization_stacks(patterns, cfg_pol);

    a = struct('id', A, 'pol', pol, 'stack1', s1, 'stack2', s2, ...
        'theta_deg', patterns(1).theta_deg, 'phi_deg', patterns(1).phi_deg, ...
        'n_el', size(s1, 1));
    % MUSIC over a 1 deg grid every step is the dominant cost. Stride 4 quarters
    % it for a measured DoA RMSE of 1.45 deg against 1.00 at stride 1 -- an
    % acceptable trade HERE because the on/off question turns on presence
    % detection rather than angular precision, and the doubled run length has to
    % be paid for somewhere. The 5 deg array is already cheap and keeps the full
    % grid.
    a.doa_stride = 1 + 3 * (mean(diff(a.theta_deg)) < 2);
    arrays.(key) = a;

    for it = 1:size(TARGETS, 1)
        th_s = TARGETS(it, 1); ph_s = TARGETS(it, 2);
        prof = kpi_array_profile(s1, s2, a.theta_deg, a.phi_deg, th_s, ph_s, ...
            config.adapt.diagonal_loading_db);

        reason = '';
        if prof.dir_s_dbi < MIN_DIR_DBI
            reason = sprintf('no gain (%.1f dBi)', prof.dir_s_dbi);
        elseif prof.guard_deg > MAX_GUARD_DEG
            reason = sprintf('beam covers sphere (guard %.0f deg)', prof.guard_deg);
        end
        fprintf('  %-22s (%3.0f,%3.0f) %5d %8.2f %8.1f %8.1f %7.3f %6d  %s\n', ...
            A, th_s, ph_s, prof.n_el, prof.dir_s_dbi, prof.hpbw_theta_deg, ...
            prof.guard_deg, prof.mirror_coh, prof.feasible_music, ...
            ternary(isempty(reason), 'test', ['SKIP: ' reason]));
        if ~isempty(reason)
            skipped{end + 1} = sprintf('%s (%.0f,%.0f): %s', A, th_s, ph_s, reason); %#ok<AGROW>
            continue
        end

        for im = 1:numel(SEP_MULT)
            sep = SEP_MULT(im) * prof.guard_deg;
            [th_j, ph_j, ok] = place_jammer(a, th_s, ph_s, sep);
            if ~ok, continue; end
            sep_true = angular_separation_deg(th_s, ph_s, th_j, ph_j);
            coh = kpi_steering_coherence(s1, s2, a.theta_deg, a.phi_deg, ...
                th_s, ph_s, th_j, ph_j);
            for ip = 1:numel(PERIODS)
                c = struct('array_id', A, 'key', key, 'target', [th_s ph_s], ...
                    'prof', prof, 'sep_mult', SEP_MULT(im), 'sep_deg', sep_true, ...
                    'theta_j_deg', th_j, 'phi_j_deg', ph_j, ...
                    'period_s', PERIODS(ip), 'coh', coh);
                cases(end + 1) = c; %#ok<AGROW>
            end
        end
    end
end

fprintf('\n%d testable cases; %d (array, target) pairs skipped as infeasible:\n', ...
    numel(cases), numel(skipped));
for i = 1:numel(skipped), fprintf('   - %s\n', skipped{i}); end

% ══════════════════════════════════════════════════════════════════
%  2. Run
% ══════════════════════════════════════════════════════════════════
timestamp  = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
for iarm = 1:numel(ARMS)
    arm = ARMS{iarm};
    out_dir = fullfile(repo_root, 'results', 'onoff_campaign', ...
        sprintf('%s_%s', timestamp, arm.tag));
    if ~isfolder(out_dir), mkdir(out_dir); end
    fprintf(['\n' repmat('=', 1, 70) '\n  ARM %d/%d: %s -> %s\n' ...
             repmat('=', 1, 70) '\n'], iarm, numel(ARMS), arm.tag, out_dir);

    csv_path = fullfile(out_dir, 'onoff_campaign.csv');
    fid = fopen(csv_path, 'w');
    if fid < 0
        error('run_onoff_campaign:CsvOpen', 'Cannot open %s.', csv_path);
    end
    fprintf(fid, ['array_id,n_el,n_comp,pol,theta_s_deg,phi_s_deg,dir_s_dbi,' ...
        'hpbw_theta_deg,guard_deg,mirror_coh,feasible_music,sep_mult,sep_deg,' ...
        'coh,theta_j_deg,phi_j_deg,period_s,duty,algorithm,' ...
        'track_score_pct,track_score_pct_std,availability_pct,oracle_gap_db,' ...
        'sinr_mean_db,oracle_sinr_mean_db,presence_err,prenull_frac,' ...
        'recovery_steps\n']);

    t_all = tic; n_fail = 0;
    for ic = 1:numel(cases)
        c = cases(ic);
        a = arrays.(c.key);
        aj = config.antijam;
        aj.theta_s_deg = c.target(1); aj.phi_s_deg = c.target(2);
        aj.sigma_s_db  = SIGMA_S_DB;  aj.jn_ratio_db = JN_RATIO_DB;
        aj.guard_deg   = c.prof.guard_deg;      % DERIVED, not the global 5 deg

        cfg_run = config;
        cfg_run.polarization = a.pol;
        cfg_run.adapt.predict.doa_stride = a.doa_stride;
        if ~isempty(arm.onoff)
            cfg_run.adapt.predict.onoff = arm.onoff;
        end

        sim_cfg = config.sim;
        sim_cfg.duration_s = CYCLES * c.period_s;

        scn_cfg = struct('id', 'ONOFF', 'motion', 'static', 'power', 'onoff', ...
            'duty_cycle', DUTY, 'toggle_period_s', c.period_s, ...
            'theta_j_deg', c.theta_j_deg, 'phi_j_deg', c.phi_j_deg, ...
            'jn_ratio_db', JN_RATIO_DB);

        acc = struct('lcmv', [], 'predict', []);
        for iseed = 1:N_SEEDS
            sim_cfg.seed = config.sim.seed + iseed - 1;
            try
                scn = sim_scenario(scn_cfg, aj, sim_cfg);
                o = closed_loop_run('oracle', a.stack1, a.stack2, a.theta_deg, ...
                    a.phi_deg, scn, aj, sim_cfg, cfg_run, []);
                for alg = {'lcmv', 'predict'}
                    g = closed_loop_run(alg{1}, a.stack1, a.stack2, a.theta_deg, ...
                        a.phi_deg, scn, aj, sim_cfg, cfg_run, []);
                    acc.(alg{1})(end + 1, :) = ...
                        onoff_metrics(g, o, scn, aj, TRACK_TOL_DB);
                end
            catch err
                n_fail = n_fail + 1;
                warning('run_onoff_campaign:CaseFailed', ...
                    '%s (%.0f,%.0f) sep %.0f period %.0f seed %d FAILED: %s (%s)', ...
                    c.array_id, c.target(1), c.target(2), c.sep_deg, c.period_s, ...
                    sim_cfg.seed, err.message, err.identifier);
            end
        end

        txt = '';
        for alg = {'lcmv', 'predict'}
            M = acc.(alg{1});
            if isempty(M), continue; end
            mu = mean(M, 1); sd = std(M, 0, 1);
            fprintf(fid, ['%s,%d,%d,%s,%.1f,%.1f,%.2f,%.1f,%.2f,%.4f,%d,' ...
                '%.2f,%.2f,%.4f,%.1f,%.1f,%.1f,%.2f,%s,' ...
                '%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.1f\n'], ...
                c.array_id, c.prof.n_el, c.prof.n_comp, a.pol, ...
                c.target(1), c.target(2), c.prof.dir_s_dbi, ...
                c.prof.hpbw_theta_deg, c.prof.guard_deg, c.prof.mirror_coh, ...
                c.prof.feasible_music, c.sep_mult, c.sep_deg, c.coh, ...
                c.theta_j_deg, c.phi_j_deg, c.period_s, DUTY, alg{1}, ...
                mu(1), sd(1), mu(2), mu(3), mu(4), mu(5), mu(6), mu(7), mu(8));
            txt = [txt sprintf(' %s %5.1f%%', alg{1}, mu(1))]; %#ok<AGROW>
        end
        el = toc(t_all);
        fprintf('  [%3d/%3d] %-22s (%3.0f,%3.0f) sep%.2fx T%.0fs %s  (%.0f s, ~%.0f s left)\n', ...
            ic, numel(cases), c.array_id, c.target(1), c.target(2), ...
            c.sep_mult, c.period_s, txt, el, el * (numel(cases) - ic) / ic);
    end
    fclose(fid);
    fprintf('\n  ARM %s complete in %.0f s (%d failed case-seeds). CSV: %s\n', ...
        arm.tag, toc(t_all), n_fail, csv_path);
end
fprintf('\nCampaign complete.\n');
end


% ────────────────────────── HELPERS ───────────────────────────────

function pol = pick_polarization(patterns)
% Use both components unless the second is numerical residue, which
% select_polarization_stacks rejects outright (it breaks MUSIC's model order).
comps = fieldnames(patterns(1).components);
c1 = stack_component(patterns, comps{1});
p1 = 10 * log10(max(abs(c1(:)).^2) + eps);
p2 = -Inf;
if numel(comps) > 1
    c2 = stack_component(patterns, comps{2});
    p2 = 10 * log10(max(abs(c2(:)).^2) + eps);
end
if (p1 - p2) > 60
    pol = comps{1};
else
    pol = 'total';
end
end


function [th_j, ph_j, ok] = place_jammer(a, th_s, ph_s, sep_deg)
% Put the jammer sep_deg from the target, preferring the theta cut and falling
% back to the phi cut when theta would run off [0, 180].
ok = true;
th_j = th_s + sep_deg;
ph_j = ph_s;
if th_j > 180 || th_j < 0
    th_j = th_s - sep_deg;
end
if th_j > 180 || th_j < 0
    % Unreachable in theta: use the phi cut, where the reachable span depends
    % on how far the target sits from a pole.
    th_j = th_s;
    if sind(th_s) < 1e-6
        ok = false; return          % at a pole, phi is degenerate
    end
    dphi = sep_deg / sind(th_s);
    if dphi > 180, ok = false; return; end
    ph_j = mod(ph_s + dphi, 360);
end
[it, ip] = nearest_index_2d(a.theta_deg, a.phi_deg, th_j, ph_j);
th_j = a.theta_deg(it);
ph_j = a.phi_deg(ip);
end


function m = onoff_metrics(g, o, scn, aj, tol_db)
% [track_score, availability, oracle_gap, sinr_mean, oracle_sinr_mean,
%  presence_err, prenull_frac, recovery_steps]
track = 100 * mean((o.sinr_db - g.sinr_db) <= tol_db);
avail = 100 * mean(g.sinr_db >= aj.sinr_min_db);
gap   = mean(o.sinr_db - g.sinr_db);

% On/off diagnostics: how well the presence indicator tracks the true duty, and
% how often the anticipatory branch actually engaged.
pres_err = NaN; prenull = NaN;
if isfield(g, 'present') && ~isempty(g.present)
    pres_err = abs(mean(g.present) - mean(scn.jammer_on));
    prenull  = mean(g.predicted_on & ~g.present);
end

% Recovery: steps from each turn-on until SINR re-crosses the threshold.
rec = recovery_steps(g.sinr_db, scn.jammer_on, aj.sinr_min_db);

m = [track, avail, gap, mean(g.sinr_db), mean(o.sinr_db), pres_err, prenull, rec];
end


function r = recovery_steps(sinr_db, on, thresh_db)
% Mean steps after a 0->1 transition before SINR is back above threshold.
turn_on = find(diff([false, on(:)']) == 1);
if isempty(turn_on), r = NaN; return; end
d = [];
for i = 1:numel(turn_on)
    k0 = turn_on(i);
    k1 = find(sinr_db(k0:end) >= thresh_db, 1, 'first');
    if isempty(k1)
        d(end + 1) = numel(sinr_db) - k0 + 1; %#ok<AGROW>
    else
        d(end + 1) = k1 - 1; %#ok<AGROW>
    end
end
r = mean(d);
end


function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
