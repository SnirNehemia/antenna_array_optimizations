function output_dir = compare_classical(config_path, n_override)
% COMPARE_CLASSICAL  Compare classical tapering/steering vs the optimizer.
%
%   output_dir = COMPARE_CLASSICAL(config_path, n_override)
%
%   MATLAB port of scripts/compare_classical.py. Reads test_config.yaml, builds
%   (synthetic) or loads (folder/CST) an N x N URA, runs each scenario through the
%   classical techniques and the optimizer, prints per-directive tables, and
%   saves a comparison figure + console log + config snapshot to a timestamped
%   results/compare_classical/ subfolder.
%
%   Inputs:
%       config_path : path to test_config.yaml. Default: <repo>/test_config.yaml.
%       n_override  : optional n_side override (N x N elements). Default: [] (use config).
%
%   Outputs:
%       output_dir : the timestamped results directory created.
%
%   Part of: Antenna Array Pattern Optimization Tool.

REPO_ROOT = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(config_path)
    config_path = fullfile(REPO_ROOT, 'test_config.yaml');
end
if nargin < 2, n_override = []; end
% Resolve a relative config name against the repo root when not in the cwd.
if ~isfile(config_path) && isfile(fullfile(REPO_ROOT, config_path))
    config_path = fullfile(REPO_ROOT, config_path);
end

cfg = load_test_config(config_path);
if ~isempty(n_override)
    cfg.n_side = n_override;
end

% ── Output directory + console log (diary) ───────────────────────
timestamp  = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
output_dir = fullfile(REPO_ROOT, 'results', 'compare_classical', timestamp);
if ~isfolder(output_dir), mkdir(output_dir); end
log_path = fullfile(output_dir, 'console_output.txt');
diary(log_path); diary on;
cleanupDiary = onCleanup(@() diary('off'));
fprintf('Output directory: %s\n', output_dir);

d_over_lambda  = cfg.d_over_lambda;
element_source = cfg.element_source;
polarization   = cfg.polarization;

% ── Element patterns ──────────────────────────────────────────────
element_patterns_secondary = [];
if strcmp(element_source, 'folder')
    [element_patterns_complex, element_patterns_secondary, theta_deg, phi_deg, n_side] = ...
        load_folder_element_patterns(cfg, REPO_ROOT);
    cfg.n_side = n_side;
    fprintf('\nLoaded %dx%d URA from folder (%d elements, %d theta x %d phi points).\n', ...
        n_side, n_side, size(element_patterns_complex, 1), numel(theta_deg), numel(phi_deg));
elseif strcmp(element_source, 'synthetic')
    if strcmp(polarization, 'total')
        error('compare_classical:TotalSynthetic', ...
            'polarization: total is not supported with element_source: synthetic.');
    end
    n_side    = cfg.n_side;
    theta_deg = (0:cfg.theta_step_deg:90)';
    phi_deg   = (0:cfg.phi_step_deg:(360 - cfg.phi_step_deg))';
    fprintf('\nBuilding %dx%d URA element patterns (%d elements, d=%g lambda, %d theta x %d phi)...\n', ...
        n_side, n_side, n_side^2, d_over_lambda, numel(theta_deg), numel(phi_deg));
    element_patterns_complex = build_ura_element_patterns(n_side, d_over_lambda, theta_deg, phi_deg);
else
    error('compare_classical:BadSource', ...
        'Unknown element_source ''%s''. Valid: synthetic, folder.', element_source);
end

% ── Scenario loop ─────────────────────────────────────────────────
n_side = cfg.n_side;
scenarios = cfg.scenarios;
all_scenario_results = cell(numel(scenarios), 1);
for si = 1:numel(scenarios)
    scenario = scenarios{si};
    fprintf('\n%s\n  %s\n%s\n', repmat('=', 1, 68), scenario_label(scenario), repmat('=', 1, 68));
    results = run_scenario(scenario, element_patterns_complex, theta_deg, phi_deg, ...
        n_side, d_over_lambda, cfg, element_patterns_secondary);
    print_metrics_table(scenario_label(scenario), results, scenario.directives);
    all_scenario_results{si} = struct('scenario', scenario, 'results', results);
end

% ── Figure ────────────────────────────────────────────────────────
figure_path = fullfile(output_dir, sprintf('compare_classical_%dx%d.png', n_side, n_side));
plot_comparison(all_scenario_results, theta_deg, n_side, cfg, figure_path);
fprintf('\nFigure saved -> %s\n', figure_path);

% ── Config snapshot ───────────────────────────────────────────────
copyfile(config_path, fullfile(output_dir, 'test_config.yaml'));
fprintf('\nFiles saved to: %s\n', output_dir);
end


% ────────────────────────── CONFIG ────────────────────────────────

function cfg = load_test_config(config_path)
% Load test config, applying built-in defaults for any missing key.
if isfile(config_path)
    cfg = read_config_yaml(config_path);
else
    fprintf('WARNING: config file not found at ''%s''. Using built-in defaults.\n', config_path);
    cfg = struct();
end
cfg = setdefault(cfg, 'n_side', 4);
cfg = setdefault(cfg, 'd_over_lambda', 0.5);
cfg = setdefault(cfg, 'theta_step_deg', 1.0);
cfg = setdefault(cfg, 'phi_step_deg', 5.0);
cfg = setdefault(cfg, 'n_restarts', 5);
cfg = setdefault(cfg, 'max_iterations', 500);
cfg = setdefault(cfg, 'cost_tolerance', 1e-9);
cfg = setdefault(cfg, 'plot_dynamic_range_db', 60.0);
cfg = setdefault(cfg, 'element_source', 'synthetic');
cfg = setdefault(cfg, 'element_patterns_dir', 'data/element_patterns/');
cfg = setdefault(cfg, 'polarization', 'copol');
if ~isfield(cfg, 'scenarios') || isempty(cfg.scenarios)
    cfg.scenarios = default_scenarios();
end
end


function s = setdefault(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end


function sc = default_scenarios()
a = struct('label', 'Scenario A - Broadside SLL suppression', ...
    'steer_theta_deg', 0.0, 'steer_phi_deg', 0.0, 'null_theta_deg', [], ...
    'directives', {{struct('type','peak','theta',7.5,'phi',180.0,'theta_width',15.0,'phi_width',360.0,'weight',1.0)}});
b = struct('label', 'Scenario B - Steered 20deg, null at 60deg', ...
    'steer_theta_deg', 20.0, 'steer_phi_deg', 0.0, 'null_theta_deg', 60.0, ...
    'directives', {{ ...
        struct('type','peak','theta',20.0,'phi',0.0,'theta_width',20.0,'phi_width',30.0,'weight',1.0), ...
        struct('type','null','theta',60.0,'phi',0.0,'theta_width',15.0,'phi_width',20.0,'weight',3.0)}});
sc = {a, b};
end


function [ep, ep_secondary, theta_deg, phi_deg, n_side] = load_folder_element_patterns(cfg, repo_root)
patterns_dir = cfg.element_patterns_dir;
if ~isfolder(patterns_dir)
    patterns_dir = fullfile(repo_root, cfg.element_patterns_dir);
end
polarization = cfg.polarization;
fprintf('  Loading element patterns from: %s\n', patterns_dir);
raw = load_element_patterns(patterns_dir);
n_elements = numel(raw);
fprintf('  Loaded %d element patterns.\n', n_elements);

n_side_f = sqrt(n_elements);
if abs(n_side_f - round(n_side_f)) > 1e-9
    error('compare_classical:NotSquare', ...
        'Number of CST files (%d) is not a perfect square.', n_elements);
end
n_side = round(n_side_f);

theta_deg = raw(1).theta_deg;
phi_deg   = raw(1).phi_deg;

ep_secondary = [];
if strcmp(polarization, 'copol')
    ep = stack_field(raw, 'E_complex');
elseif strcmp(polarization, 'cross')
    ep = stack_field(raw, 'cross_complex');
elseif strcmp(polarization, 'total')
    ep           = stack_field(raw, 'E_complex');
    ep_secondary = stack_field(raw, 'cross_complex');
    peak_amp = max([max(abs(ep), [], 'all'), max(abs(ep_secondary), [], 'all')]);
    if peak_amp > 1e-30
        ep = ep / peak_amp;
        ep_secondary = ep_secondary / peak_amp;
    end
    return;
else
    error('compare_classical:BadPolarization', ...
        'Unknown polarization ''%s''. Valid: copol, cross, total.', polarization);
end

peak_amp = max(abs(ep), [], 'all');
if peak_amp > 1e-30
    ep = ep / peak_amp;
end
end


function stack = stack_field(patterns, key)
n = numel(patterns);
[nt, np_] = size(patterns(1).(key));
stack = zeros(n, nt, np_);
for i = 1:n
    stack(i, :, :) = patterns(i).(key);
end
end


% ────────────────────────── CONSOLE TABLE ─────────────────────────

function print_metrics_table(scenario_label_str, results, directives)
CLASSICAL = {'uniform','hamming','hanning','kaiser_b3','kaiser_b6', ...
    'chebyshev_25dB','chebyshev_40dB','taylor_25dB','optimized'};
techs = CLASSICAL(cellfun(@(nm) isfield(results, nm), CLASSICAL));

fprintf('\n%s\n  %s\n%s\n', repmat('=', 1, 75), scenario_label_str, repmat('=', 1, 75));
for di = 1:numel(directives)
    d = directives{di};
    sym = get_field(d, 'width', NaN);
    fprintf('\n  Directive %d: %s  theta=%.0f phi=%.0f theta_width=%s\n', ...
        di, upper(d.type), d.theta, get_field(d, 'phi', 0.0), ...
        num2str(get_field(d, 'theta_width', sym)));
    fprintf('  %-18s %11s %10s %10s %10s\n', 'Technique', 'Peak(dBi)', 'Min(dB)', 'Mean(dB)', 'Max(dB)');
    for ti = 1:numel(techs)
        dm = results.(techs{ti}).directive_metrics(di);
        marker = '';
        if strcmp(techs{ti}, 'optimized'), marker = ' <--'; end
        fprintf('  %-18s %11.1f %10.1f %10.1f %10.1f%s\n', ...
            techs{ti}, dm.peak_dbi, dm.min_db, dm.mean_db, dm.max_db, marker);
    end
end
end


function s = scenario_label(scenario)
if isfield(scenario, 'label') && ~isempty(scenario.label)
    s = scenario.label;
else
    s = 'scenario';
end
end


function val = get_field(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end
