function output_dir = run_jammer_demo(config_path)
% RUN_JAMMER_DEMO  Scenario demo runner: closed loops + KPI table + GIFs.
%
%   output_dir = RUN_JAMMER_DEMO(config_path)
%
%   Driven entirely by jammer_config.yaml (same section structure as
%   config.yaml plus a `gif` section): the user dictates the jammer
%   scenarios there — static on/off at fixed angles/amplitudes, power steps,
%   or constant-speed motion (sim_scenario vocabulary incl. the angle_deg
%   and per-scenario jn_ratio_db / duration_s overrides).
%
%   For every scenario x algorithm in antijam.algorithms (single seed —
%   this is a demo runner, not the Monte Carlo campaign):
%     - runs the closed loop (closed_loop_run) with the oracle timeline as
%       reference,
%     - evaluates the 5 KPIs (kpi_evaluate -> kpi_table.csv/txt),
%     - renders vid_<scenario>_<algorithm>.<ext> (save_run_gif): the full 2-D
%       radiation pattern evolving over time, the jammer's position/strength
%       as a dot, and live SINR + gain-to-target traces. Extension follows
%       gif.format ('gif' default, or 'mp4' for MPEG-4).
%
%   Inputs:
%       config_path : path to the demo YAML. Default: <repo>/jammer_config.yaml.
%
%   Outputs:
%       output_dir : the timestamped results/jammer_demo/<ts>/ folder.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

REPO_ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1 || isempty(config_path)
    config_path = fullfile(REPO_ROOT, 'jammer_config.yaml');
end
if ~isfile(config_path) && isfile(fullfile(REPO_ROOT, config_path))
    config_path = fullfile(REPO_ROOT, config_path);
end

% ── Config ────────────────────────────────────────────────────────
config = read_config_yaml(config_path);
for key = {'element_patterns_dir', 'polarization', 'antijam', 'sim', 'adapt', ...
           'agent', 'gif', 'output'}
    if ~isfield(config, key{1})
        error('run_jammer_demo:MissingKey', 'Missing required config key: ''%s''.', key{1});
    end
end
aj = config.antijam;
for key = {'theta_s_deg', 'phi_s_deg', 'sinr_min_db', 'algorithms', 'scenarios'}
    if ~isfield(aj, key{1}) || isempty(aj.(key{1}))
        error('run_jammer_demo:MissingKey', ...
            'Missing required antijam config key: ''%s''.', key{1});
    end
end
algorithms = aj.algorithms;

% ── Output folder ─────────────────────────────────────────────────
results_dir = config.output.results_dir;
if ~(ispc && ~isempty(regexp(results_dir, '^([A-Za-z]:|\\\\)', 'once'))) && results_dir(1) ~= '/'
    results_dir = fullfile(REPO_ROOT, results_dir);
end
timestamp  = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
output_dir = fullfile(results_dir, 'jammer_demo', timestamp);
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

% ── Codebook (cached, shared with run_antijam) ────────────────────
codebook = [];
if any(strcmp(algorithms, 'bandit'))
    if ~isfield(config.agent, 'optimizer')
        error('run_jammer_demo:MissingKey', ...
            'Missing required agent config key: ''optimizer'' (codebook synthesis settings).');
    end
    cache_path = fullfile(results_dir, 'antijam', 'codebook_cache.mat');
    codebook = agent_codebook_build(stack1, stack2, theta_deg, phi_deg, ...
        aj, config.agent, config.agent.optimizer, cache_path);
end

% ── Scenarios x algorithms: run, evaluate, render ─────────────────
runs = {};
for isc = 1:numel(aj.scenarios)
    scn_cfg = aj.scenarios{isc};
    scn = sim_scenario(scn_cfg, aj, config.sim);
    fprintf('%s (theta_j start %.0f deg): ', scn.id, scn.theta_j_deg(1));

    o_log = closed_loop_run('oracle', stack1, stack2, theta_deg, phi_deg, scn, aj, config.sim, ...
                            config, codebook);
    o_log.oracle_sinr_db = o_log.sinr_db;

    for ia = 1:numel(algorithms)
        alg = algorithms{ia};
        if strcmp(alg, 'oracle')
            log = o_log;
        else
            log = closed_loop_run(alg, stack1, stack2, theta_deg, phi_deg, scn, aj, config.sim, ...
                                  config, codebook);
            log.oracle_sinr_db = o_log.sinr_db;
        end
        kpi = kpi_evaluate(log, scn, aj);
        runs{end + 1} = struct('scenario_id', scn.id, 'algorithm', alg, ...
            'seed', config.sim.seed, 'run_log', log, 'kpi', kpi, ...
            'scenario', scn); %#ok<AGROW>
        vid_ext = 'gif';
        if isfield(config.gif, 'format') && ~isempty(config.gif.format)
            vid_ext = lower(config.gif.format);
        end
        vid_path = fullfile(output_dir, sprintf('vid_%s_%s.%s', scn.id, alg, vid_ext));
        save_run_gif(log, scn, stack1, stack2, theta_deg, phi_deg, aj, ...
            config.gif, sprintf('%s / %s', scn.id, alg), vid_path);
        fprintf('%s ', alg);
    end
    fprintf('\n');
end

write_kpi_table(runs, aj, output_dir);
fprintf('Demo written to %s\n', output_dir);
end
