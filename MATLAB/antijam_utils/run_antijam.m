function output_dir = run_antijam(config_path)
% RUN_ANTIJAM  Closed-loop anti-jam campaign: scenarios x algorithms -> report.
%
%   output_dir = RUN_ANTIJAM(config_path)
%
%   Top-level driver (analogue of run_optimization for this milestone):
%   loads config + CST element patterns, extracts the 1-D azimuth cut,
%   builds/loads the codebook (cached), runs every scenario x algorithm x
%   Monte Carlo seed closed loop, evaluates the 5 KPIs per run, and writes a
%   KPI table (CSV + txt) + report figures + config snapshot to a timestamped
%   results/antijam/<ts>/ folder.
%
%   Modes are fixed per algorithm: oracle, lcmv -> Mode C; spsa, bandit ->
%   Mode S. The oracle runs once per scenario x seed and its SINR timeline is
%   shared as the oracle reference of every other algorithm's log (the scalar
%   SINR channel is deterministic given the scenario, so timelines align).
%   The oracle applies the perfect-knowledge LCMV weights with a one-step lag
%   (R of step k -> w applied at k+1) — negligible at the configured drift.
%
%   Inputs:
%       config_path : path to config.yaml. Default: <repo>/config.yaml.
%
%   Outputs:
%       output_dir : the timestamped results directory that was created.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

REPO_ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1 || isempty(config_path)
    config_path = fullfile(REPO_ROOT, 'config.yaml');
end
if ~isfile(config_path) && isfile(fullfile(REPO_ROOT, config_path))
    config_path = fullfile(REPO_ROOT, config_path);
end

% ── Config ────────────────────────────────────────────────────────
config = read_config_yaml(config_path);
for key = {'element_patterns_dir', 'antijam', 'sim', 'adapt', 'agent', 'output'}
    if ~isfield(config, key{1})
        error('run_antijam:MissingKey', 'Missing required config key: ''%s''.', key{1});
    end
end
aj  = config.antijam;
smc = config.sim;
for key = {'theta_s_deg', 'sinr_min_db', 'algorithms', 'scenarios', 'cut_type'}
    if ~isfield(aj, key{1}) || isempty(aj.(key{1}))
        error('run_antijam:MissingKey', 'Missing required antijam config key: ''%s''.', key{1});
    end
end
if ~isfield(smc, 'n_runs') || isempty(smc.n_runs)
    error('run_antijam:MissingKey', 'Missing required sim config key: ''n_runs''.');
end
algorithms = aj.algorithms;

% ── Output folder ─────────────────────────────────────────────────
results_dir = config.output.results_dir;
if ~(ispc && ~isempty(regexp(results_dir, '^([A-Za-z]:|\\\\)', 'once'))) && results_dir(1) ~= '/'
    results_dir = fullfile(REPO_ROOT, results_dir);
end
timestamp  = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
output_dir = fullfile(results_dir, 'antijam', timestamp);
if ~isfolder(output_dir), mkdir(output_dir); end
copyfile(config_path, fullfile(output_dir, 'config_snapshot.yaml'));
fprintf('Output directory: %s\n', output_dir);

% ── Element patterns + polarization (Milestone-1 convention) ──────
patterns_dir = config.element_patterns_dir;
if ~isfolder(patterns_dir)
    patterns_dir = fullfile(REPO_ROOT, patterns_dir);
end
fprintf('Loading element patterns...\n');
patterns  = load_element_patterns(patterns_dir);
theta_deg = patterns(1).theta_deg;
phi_deg   = patterns(1).phi_deg;
[stack1, stack2, pol_label] = select_polarization_stacks(patterns, config);
fprintf('  %d elements, grid %dx%d, polarization: %s\n', ...
    size(stack1, 1), numel(theta_deg), numel(phi_deg), pol_label);

% ── 1-D azimuth cut ───────────────────────────────────────────────
[E1c, E2c, cut_ang] = extract_cut(stack1, stack2, theta_deg, phi_deg, aj);
fprintf('  Cut: %s, %d angles\n', aj.cut_type, numel(cut_ang));

% ── Codebook (cached) ─────────────────────────────────────────────
codebook = [];
if any(strcmp(algorithms, 'bandit'))
    if ~isfield(config.agent, 'optimizer')
        error('run_antijam:MissingKey', ...
            'Missing required agent config key: ''optimizer'' (codebook synthesis settings).');
    end
    cache_path = fullfile(results_dir, 'antijam', 'codebook_cache.mat');
    codebook = agent_codebook_build(stack1, stack2, theta_deg, phi_deg, ...
        aj, config.agent, config.agent.optimizer, cache_path);
end

% ── Campaign: scenarios x seeds x algorithms ──────────────────────
runs = {};
for isc = 1:numel(aj.scenarios)
    scn_cfg = aj.scenarios{isc};
    for run_i = 1:smc.n_runs
        sim_cfg = smc;
        sim_cfg.seed = smc.seed + run_i - 1;
        scn = sim_scenario(scn_cfg, aj, sim_cfg);
        fprintf('%s seed %d: ', scn.id, sim_cfg.seed);

        % Oracle first — its SINR timeline is the reference for all logs.
        o_log = closed_loop('oracle', E1c, E2c, cut_ang, scn, aj, sim_cfg, config, codebook);
        o_log.oracle_sinr_db = o_log.sinr_db;

        for ia = 1:numel(algorithms)
            alg = algorithms{ia};
            if strcmp(alg, 'oracle')
                log = o_log;
            else
                log = closed_loop(alg, E1c, E2c, cut_ang, scn, aj, sim_cfg, config, codebook);
                log.oracle_sinr_db = o_log.sinr_db;
            end
            kpi = kpi_evaluate(log, scn, aj);
            runs{end + 1} = struct('scenario_id', scn.id, 'algorithm', alg, ...
                'seed', sim_cfg.seed, 'run_log', log, 'kpi', kpi, ...
                'scenario', scn); %#ok<AGROW>
            fprintf('%s ', alg);
        end
        fprintf('\n');
    end
end

% ── KPI table + report ────────────────────────────────────────────
write_kpi_table(runs, aj, output_dir);
plot_antijam_report(runs, config, output_dir);
if ~isempty(codebook)
    save(fullfile(output_dir, 'codebook.mat'), 'codebook');
end
fprintf('Report written to %s\n', output_dir);
end


% ────────────────────────── HELPERS ───────────────────────────────

function [stack1, stack2, label] = select_polarization_stacks(patterns, config)
% Mirror run_optimization's polarization convention: named component or
% 'total' (incoherent sum of exactly two detected components).
if ~isfield(config, 'polarization') || isempty(config.polarization)
    error('run_antijam:MissingKey', 'Missing required config key: ''polarization''.');
end
pol   = config.polarization;
names = sort(fieldnames(patterns(1).components));
if strcmpi(pol, 'total')
    if numel(names) ~= 2
        error('run_antijam:BadPolarization', ...
            'polarization ''total'' requires exactly 2 components, found {%s}.', ...
            strjoin(names, ', '));
    end
    stack1 = stack_component(patterns, names{1});
    stack2 = stack_component(patterns, names{2});
    label  = sprintf('total (%s + %s)', names{1}, names{2});
else
    hit = names(strcmpi(names, pol));
    if isempty(hit)
        error('run_antijam:BadPolarization', ...
            'polarization ''%s'' not found among components {%s}.', ...
            pol, strjoin(names, ', '));
    end
    stack1 = stack_component(patterns, hit{1});
    stack2 = [];
    label  = hit{1};
end
end


function [E1c, E2c, cut_ang] = extract_cut(stack1, stack2, theta_deg, phi_deg, aj)
% Extract the engine's 1-D azimuth cut per antijam.cut_type.
switch aj.cut_type
    case 'phi_cut'
        if ~isfield(aj, 'cut_theta_deg') || isempty(aj.cut_theta_deg)
            error('run_antijam:MissingKey', 'phi_cut requires antijam.cut_theta_deg.');
        end
        it = nearest_index(theta_deg(:), aj.cut_theta_deg);
        cut_ang = phi_deg(:).';
        E1c = squeeze(stack1(:, it, :));
        E2c = [];
        if ~isempty(stack2), E2c = squeeze(stack2(:, it, :)); end
    case 'theta_cut'
        if ~isfield(aj, 'cut_phi_deg') || isempty(aj.cut_phi_deg)
            error('run_antijam:MissingKey', 'theta_cut requires antijam.cut_phi_deg.');
        end
        ip0  = nearest_index(phi_deg(:), mod(aj.cut_phi_deg, 360));
        ip180 = nearest_index(phi_deg(:), mod(aj.cut_phi_deg + 180, 360));
        cut_ang = [-fliplr(theta_deg(2:end)), theta_deg];
        E1c = [fliplr(squeeze(stack1(:, 2:end, ip180))), squeeze(stack1(:, :, ip0))];
        E2c = [];
        if ~isempty(stack2)
            E2c = [fliplr(squeeze(stack2(:, 2:end, ip180))), squeeze(stack2(:, :, ip0))];
        end
    otherwise
        error('run_antijam:BadCutType', 'Unknown cut_type ''%s''.', aj.cut_type);
end
end


function log = closed_loop(alg, E1c, E2c, cut_ang, scn, aj, sim_cfg, config, codebook)
% One closed-loop run of `alg` against `scn`. Returns the run log.
mode = 'S';
if any(strcmp(alg, {'oracle', 'lcmv'})), mode = 'C'; end
st   = sim_engine_init(E1c, E2c, cut_ang, scn, aj, sim_cfg, mode);
n_el = size(E1c, 1);
T    = numel(scn.t_s);

switch alg
    case 'oracle'
        w = adapt_lcmv(eye(n_el), st.e_s, 0);           % quiescent start
    case 'lcmv'
        trk = adapt_tracking_init(config.adapt, st.e_s, n_el);
        w = trk.w;
    case 'spsa'
        % Warm start from the quiescent MVDR beam — the uniform start sits
        % ~10-20 dB deeper on the real array and dominates the probe budget.
        w_q = adapt_lcmv(eye(n_el), st.e_s, 0);
        sp = adapt_spsa_init(config.adapt.spsa, w_q, sim_cfg.seed + 5000);
        w = sp.w;
    case 'bandit'
        bd = agent_bandit_init(config.agent, codebook, sim_cfg.seed + 6000);
        w = bd.w;
    otherwise
        error('run_antijam:BadAlgorithm', 'Unknown algorithm ''%s''.', alg);
end

log = struct();
log.sinr_db   = zeros(1, T);
log.W         = complex(zeros(n_el, T));
log.arm_index = NaN(1, T);
for k = 1:T
    log.W(:, k) = w;
    [obs, st] = sim_engine_step(st, w);
    log.sinr_db(k) = obs.sinr_db;
    switch alg
        case 'oracle'
            w = adapt_lcmv(sim_analytic_covariance(st), st.e_s, 0);
        case 'lcmv'
            [w, trk] = adapt_tracking_update(trk, obs);
        case 'spsa'
            [w, sp] = adapt_spsa_update(sp, obs);
        case 'bandit'
            log.arm_index(k) = bd.arm_index;
            [w, bd] = agent_bandit_update(bd, obs);
    end
end
log.cut = struct('E1', E1c, 'E2', E2c, 'angle_deg', cut_ang, 'e_s', st.e_s);
if strcmp(alg, 'bandit')
    log.null_center_deg = codebook.null_center_deg;   % for the arm heatmap
end
end


function write_kpi_table(runs, aj, output_dir)
% Aggregate per scenario x algorithm; write CSV + printable txt table.
key = cellfun(@(r) sprintf('%s|%s', r.scenario_id, r.algorithm), runs, ...
              'UniformOutput', false);
[groups, ~, gidx] = unique(key, 'stable');
lines_csv = {['scenario,algorithm,n_runs,availability_mean,recovery_median_steps,' ...
              'null_err_median_deg,gain_penalty_mean_db,oracle_gap_mean_db']};
lines_txt = {sprintf('%-4s %-8s %6s %10s %10s %10s %10s %10s', ...
    'scn', 'alg', 'runs', 'avail', 'recov', 'null_err', 'gain_pen', 'ora_gap')};
for g = 1:numel(groups)
    members = runs(gidx == g);
    parts = strsplit(groups{g}, '|');
    avail = mean(cellfun(@(r) r.kpi.availability, members));
    recs  = cellfun(@(r) {r.kpi.recovery_steps}, members);
    recs  = [recs{:}];
    rec_med = median(recs(~isnan(recs)));
    if isempty(rec_med) || isnan(rec_med), rec_med = NaN; end
    nerr  = cellfun(@(r) {r.kpi.null_pointing_err_deg}, members);
    nerr  = [nerr{:}];
    nerr_med = median(nerr(~isnan(nerr)));
    gpen  = mean(cellfun(@(r) mean(r.kpi.peak_gain_penalty_db), members));
    ogap  = mean(cellfun(@(r) r.kpi.oracle_gap_mean_db, members));
    lines_csv{end + 1} = sprintf('%s,%s,%d,%.4f,%.1f,%.2f,%.2f,%.2f', ...
        parts{1}, parts{2}, numel(members), avail, rec_med, nerr_med, gpen, ogap); %#ok<AGROW>
    lines_txt{end + 1} = sprintf('%-4s %-8s %6d %10.3f %10.1f %10.2f %10.2f %10.2f', ...
        parts{1}, parts{2}, numel(members), avail, rec_med, nerr_med, gpen, ogap); %#ok<AGROW>
end
fid = fopen(fullfile(output_dir, 'kpi_table.csv'), 'w');
fprintf(fid, '%s\n', lines_csv{:});
fclose(fid);
fid = fopen(fullfile(output_dir, 'kpi_table.txt'), 'w');
fprintf(fid, 'Anti-jam campaign KPI table (threshold %.1f dB)\n\n', aj.sinr_min_db);
fprintf(fid, '%s\n', lines_txt{:});
fclose(fid);
fprintf('%s\n', lines_txt{:});
end
