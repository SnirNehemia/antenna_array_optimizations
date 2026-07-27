function output_dir = run_antijam(config_path)
% RUN_ANTIJAM  Closed-loop anti-jam campaign: scenarios x algorithms -> report.
%
%   output_dir = RUN_ANTIJAM(config_path)
%
%   Top-level driver (analogue of run_optimization for this milestone):
%   loads config + CST element patterns, [P7] operates on the full 2-D
%   (theta, phi) grid natively (no 1-D cut extraction), builds/loads the
%   codebook (cached), runs every scenario x algorithm x
%   Monte Carlo seed closed loop (closed_loop_run), evaluates the 5 KPIs per
%   run, and writes a KPI table (CSV + txt) + report figures + config
%   snapshot to a timestamped results/antijam/<ts>/ folder.
%
%   The oracle runs once per scenario x seed and its SINR timeline is shared
%   as the oracle reference of every other algorithm's log (the scalar SINR
%   channel is deterministic given the scenario, so timelines align).
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
for key = {'theta_s_deg', 'phi_s_deg', 'sinr_min_db', 'algorithms', 'scenarios'}
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

% ── Element patterns, polarization, cut ───────────────────────────
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
        o_log = closed_loop_run('oracle', stack1, stack2, theta_deg, phi_deg, scn, aj, sim_cfg, config, codebook);
        o_log.oracle_sinr_db = o_log.sinr_db;

        for ia = 1:numel(algorithms)
            alg = algorithms{ia};
            if strcmp(alg, 'oracle')
                log = o_log;
            else
                log = closed_loop_run(alg, stack1, stack2, theta_deg, phi_deg, scn, aj, sim_cfg, config, codebook);
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
