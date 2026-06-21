function output_dir = run_optimization(config_path)
% RUN_OPTIMIZATION  Top-level pipeline: load config, optimize, save results.
%
%   output_dir = RUN_OPTIMIZATION(config_path)
%
%   MATLAB port of scripts/run_optimization.py. Loads a YAML config, parses CST
%   element patterns, runs the multi-start fmincon optimizer, evaluates metrics,
%   saves plots (+ optional GIF), and writes weights.csv / metrics.json /
%   run_report.txt / a config snapshot to a timestamped results subfolder.
%
%   Inputs:
%       config_path : path to config.yaml. Default: <repo>/config.yaml.
%
%   Outputs:
%       output_dir : the timestamped results directory that was created.
%
%   Part of: Antenna Array Pattern Optimization Tool.

REPO_ROOT = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(config_path)
    config_path = fullfile(REPO_ROOT, 'config.yaml');
end
% A relative config name not found in the current folder is resolved against the
% repo root, so the script works regardless of the working directory.
if ~isfile(config_path) && isfile(fullfile(REPO_ROOT, config_path))
    config_path = fullfile(REPO_ROOT, config_path);
end

% Special polarization keyword: incoherent power sum across all detected
% polarization components (e.g. |AF_Copol|^2 + |AF_Cross|^2, or
% |AF_Theta|^2 + |AF_Phi|^2). Only supported when exactly two components are
% present, since the optimizer/cost function take a single secondary stack.
POLARIZATION_TOTAL = 'total';

% Default polarization when not specified in config.yaml.
DEFAULT_POLARIZATION = 'copol';

REQUIRED_CONFIG_KEYS = {'element_patterns_dir', 'directives', 'optimizer', 'output'};

% ── Load and validate config ──────────────────────────────────────
config = read_config_yaml(config_path);
for i = 1:numel(REQUIRED_CONFIG_KEYS)
    if ~isfield(config, REQUIRED_CONFIG_KEYS{i})
        error('run_optimization:MissingKey', ...
            'Missing required config key: ''%s''.', REQUIRED_CONFIG_KEYS{i});
    end
end
polarization = get_field(config, 'polarization', DEFAULT_POLARIZATION);

% ── Timestamped output directory ──────────────────────────────────
results_dir = resolve_path(config.output.results_dir, REPO_ROOT);
timestamp   = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
output_dir  = fullfile(results_dir, 'optimizations', timestamp);
if ~isfolder(output_dir), mkdir(output_dir); end
fprintf('Output directory: %s\n', output_dir);

% ── Stage 1: load element patterns ────────────────────────────────
fprintf('Loading element patterns...\n');
patterns_dir = resolve_path(config.element_patterns_dir, REPO_ROOT);
patterns  = load_element_patterns(patterns_dir);
theta_deg = patterns(1).theta_deg;
phi_deg   = patterns(1).phi_deg;
fprintf('  Polarization: %s\n', polarization);

% Component names are matched case-insensitively against whatever
% Abs(<name>)/Phase(<name>) pairs were detected in the CST export header
% (e.g. 'Copol'/'Cross' or 'Theta'/'Phi').
component_names = sort(fieldnames(patterns(1).components));

if strcmpi(polarization, POLARIZATION_TOTAL)
    % Total field: incoherent power sum across all detected components,
    % e.g. |AF_Copol|^2 + |AF_Cross|^2 or |AF_Theta|^2 + |AF_Phi|^2.
    if numel(component_names) ~= 2
        error('run_optimization:BadPolarization', ...
            'polarization: ''total'' requires exactly 2 polarization components, found {%s}.', ...
            strjoin(component_names, ', '));
    end
    name_primary   = component_names{1};
    name_secondary = component_names{2};
    element_patterns_stacked   = stack_component(patterns, name_primary);
    element_patterns_secondary = stack_component(patterns, name_secondary);
    fprintf('  Using ''%s'' + ''%s'' -> total power sum.\n', name_primary, name_secondary);
else
    matches = component_names(strcmpi(component_names, polarization));
    if isempty(matches)
        error('run_optimization:BadPolarization', ...
            'Polarization ''%s'' not found. Available components: {%s} (or ''total'').', ...
            polarization, strjoin(component_names, ', '));
    end
    fprintf('  Using component ''%s''\n', matches{1});
    element_patterns_stacked   = stack_component(patterns, matches{1});
    element_patterns_secondary = [];
end
[n_elements, n_theta, n_phi] = size(element_patterns_stacked);
fprintf('  Loaded %d elements, grid %dx%d\n', n_elements, n_theta, n_phi);

% ── Stage 4: optimizer ────────────────────────────────────────────
fprintf('Running optimizer...\n');
directives = config.directives;
t_start    = tic;
opt = run_optimizer(element_patterns_stacked, theta_deg, phi_deg, ...
    directives, config.optimizer, element_patterns_secondary);
elapsed_sec = toc(t_start);

weights_complex = opt.weights_complex;
cost_history    = opt.cost_history;

% ── Power grid (metrics + plotting) ───────────────────────────────
af_primary = compute_array_factor(weights_complex, element_patterns_stacked);
if ~isempty(element_patterns_secondary)
    af_secondary = compute_array_factor(weights_complex, element_patterns_secondary);
    power_linear = abs(af_primary) .^ 2 + abs(af_secondary) .^ 2;
    metrics_af   = sqrt(power_linear);   % synthetic real AF: |sqrt(P)|^2 == P
else
    power_linear = abs(af_primary) .^ 2;
    metrics_af   = [];
end
array_factor_db_grid = 10 * log10(max(power_linear, 1e-30));

% ── Stage 5: metrics ──────────────────────────────────────────────
metrics = evaluate_metrics(element_patterns_stacked, theta_deg, phi_deg, ...
    weights_complex, directives, cost_history, metrics_af);

% ── Stage 6: plots ────────────────────────────────────────────────
fprintf('Saving plots...\n');
save_all_plots(array_factor_db_grid, theta_deg, phi_deg, weights_complex, directives, ...
    opt.all_cost_histories, opt.best_run_index, opt.all_run_labels, config.output, output_dir);

% ── Stage 7: pattern GIF (optional) ──────────────────────────────
if get_field(config.output, 'save_pattern_gif', false)
    fprintf('Saving pattern evolution GIF...\n');
    save_pattern_gif(opt.weights_history, element_patterns_stacked, theta_deg, phi_deg, ...
        directives, cost_history, config.output, output_dir, element_patterns_secondary);
end

% ── Save weights / metrics / report / config snapshot ────────────
save_weights_csv(weights_complex, fullfile(output_dir, 'weights.csv'));
save_metrics_json(metrics, fullfile(output_dir, 'metrics.json'));
save_run_report(metrics, opt, elapsed_sec, config, fullfile(output_dir, 'run_report.txt'));
copyfile(config_path, fullfile(output_dir, 'config.yaml'));

% ── Summary ───────────────────────────────────────────────────────
fprintf('\nResults:\n');
fprintf('  Optimization time:  %.1f s\n', elapsed_sec);
fprintf('  Converged:          %d\n', opt.result.success);
fprintf('  Iterations:         %d\n', metrics.iteration_count);
fprintf('  Final cost J:       %.6f\n', metrics.total_cost);
fprintf('  Global peak:        %.2f dBi\n', metrics.global_peak_dbi);
fprintf('\nFiles saved to: %s\n', output_dir);
end


% ────────────────────────── HELPERS ───────────────────────────────

function val = get_field(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end


function p = resolve_path(p, repo_root)
% Resolve a possibly-relative config path against the repo root. Absolute or
% already-existing paths are returned unchanged (do NOT prepend repo_root).
if isfolder(p) || isfile(p) || is_absolute_path(p)
    return;
end
p = fullfile(repo_root, p);
end


function tf = is_absolute_path(p)
% True for a Windows drive path (C:\...), UNC (\\...), or POSIX-rooted (/...).
tf = ~isempty(regexp(p, '^([A-Za-z]:[\\/]|[\\/])', 'once'));
end


function save_weights_csv(weights_complex, output_path)
% Columns: element_index, amplitude, phase_deg, real, imag (2-decimal like Python).
fid = fopen(output_path, 'w');
fprintf(fid, 'element_index,amplitude,phase_deg,real,imag\n');
for n = 1:numel(weights_complex)
    w = weights_complex(n);
    fprintf(fid, '%d,%.2f,%.2f,%.2f,%.2f\n', ...
        n - 1, abs(w), angle(w) * 180 / pi, real(w), imag(w));
end
fclose(fid);
end


function save_metrics_json(metrics, output_path)
% Convert empty ([] = Python None) to NaN so jsonencode emits null.
m = metrics;
for k = 1:numel(m.directive_metrics)
    if isempty(m.directive_metrics(k).null_depth_db)
        m.directive_metrics(k).null_depth_db = NaN;
    end
end
if isempty(m.peak_to_null_ratio_db)
    m.peak_to_null_ratio_db = NaN;
end
fid = fopen(output_path, 'w');
fprintf(fid, '%s', jsonencode(m, 'PrettyPrint', true));
fclose(fid);
end


function save_run_report(metrics, opt, elapsed_sec, config, output_path)
opt_cfg    = config.optimizer;
directives = config.directives;
lines = {};
bar = repmat('=', 1, 56);
lines{end+1} = bar;
lines{end+1} = 'Antenna Array Optimization Report';
lines{end+1} = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
lines{end+1} = bar;
lines{end+1} = '';
lines{end+1} = 'RUN INFO';
lines{end+1} = sprintf('  Elapsed time:      %.1f s', elapsed_sec);
lines{end+1} = sprintf('  Converged:         %d', opt.result.success);
lines{end+1} = sprintf('  Iterations:        %d', metrics.iteration_count);
lines{end+1} = sprintf('  Final cost J:      %.6f', metrics.total_cost);
lines{end+1} = '';
lines{end+1} = 'OPTIMIZER SETTINGS';
lines{end+1} = '  Method:            fmincon (interior-point)';
lines{end+1} = sprintf('  Max iterations:    %s', num2str(get_field(opt_cfg, 'max_iterations', NaN)));
lines{end+1} = sprintf('  Cost tolerance:    %g', get_field(opt_cfg, 'cost_tolerance', NaN));
lines{end+1} = sprintf('  N restarts:        %s', num2str(get_field(opt_cfg, 'n_restarts', NaN)));
lines{end+1} = '';
lines{end+1} = 'DIRECTIVES';
for i = 1:numel(directives)
    d = directives{i};
    sym = get_field(d, 'width', NaN);
    lines{end+1} = sprintf('  [%d] %-4s  theta=%.1f  phi=%.1f  theta_width=%s  phi_width=%s  weight=%g', ...
        i - 1, d.type, d.theta, get_field(d, 'phi', 0.0), ...
        num2str(get_field(d, 'theta_width', sym)), num2str(get_field(d, 'phi_width', sym)), ...
        get_field(d, 'weight', 1.0));
end
lines{end+1} = '';
lines{end+1} = 'RESULTS';
lines{end+1} = sprintf('  Global peak:       %.2f dBi', metrics.global_peak_dbi);
if ~isempty(metrics.peak_to_null_ratio_db) && ~isnan(metrics.peak_to_null_ratio_db)
    lines{end+1} = sprintf('  Peak-to-null:      %.2f dB', metrics.peak_to_null_ratio_db);
end
lines{end+1} = '';
for i = 1:numel(metrics.directive_metrics)
    dm = metrics.directive_metrics(i);
    if strcmp(dm.type, 'peak')
        lines{end+1} = sprintf('  Peak @ theta=%.1f, phi=%.1f:  %.2f dBi', ...
            dm.theta_deg, dm.phi_deg, dm.gain_dbi);
    else
        lines{end+1} = sprintf('  Null @ theta=%.1f, phi=%.1f:  %.2f dBi  (depth = %.2f dB)', ...
            dm.theta_deg, dm.phi_deg, dm.gain_dbi, dm.null_depth_db);
    end
end
lines{end+1} = '';
lines{end+1} = bar;

fid = fopen(output_path, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
