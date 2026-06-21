function metrics = manual_weights_render(config_path, weights_csv_path, output_path)
% MANUAL_WEIGHTS_RENDER  Non-interactive replacement for the tkinter weight tuner.
%
%   metrics = MANUAL_WEIGHTS_RENDER(config_path, weights_csv_path, output_path)
%
%   Loads a config + a weights CSV (the format written by run_optimization), computes
%   the array factor / directivity / metrics for the configured directives, saves a
%   2-D directivity (dBi) heatmap with directive overlays, and prints a metric
%   summary. Replaces scripts/manual_weights.py's live GUI with a batch render
%   (the interactive GUI is out of scope per the agreed plan).
%
%   Inputs:
%       config_path      : path to config.yaml (needs element_patterns_dir, directives).
%       weights_csv_path : optional weights CSV (amplitude, phase_deg columns).
%                          [] or omitted -> uniform weights.
%       output_path      : optional PNG path. Default: <repo>/results/manual_render.png.
%
%   Outputs:
%       metrics : the evaluate_metrics struct for the rendered pattern.
%
%   Part of: Antenna Array Pattern Optimization Tool.

REPO_ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1 || isempty(config_path)
    config_path = fullfile(REPO_ROOT, 'config.yaml');
end
% Resolve a relative config name against the repo root when not in the cwd.
if ~isfile(config_path) && isfile(fullfile(REPO_ROOT, config_path))
    config_path = fullfile(REPO_ROOT, config_path);
end
if nargin < 2, weights_csv_path = []; end
if nargin < 3 || isempty(output_path)
    output_path = fullfile(REPO_ROOT, 'results', 'manual_render.png');
end
out_parent = fileparts(output_path);
if ~isempty(out_parent) && ~isfolder(out_parent), mkdir(out_parent); end

config = read_config_yaml(config_path);
for key = {'element_patterns_dir', 'directives'}
    if ~isfield(config, key{1})
        error('manual_weights_render:MissingKey', 'Required config key ''%s'' missing.', key{1});
    end
end

patterns_dir = config.element_patterns_dir;
if ~isfolder(patterns_dir)
    patterns_dir = fullfile(REPO_ROOT, config.element_patterns_dir);
end
patterns  = load_element_patterns(patterns_dir);
theta_deg = patterns(1).theta_deg;
phi_deg   = patterns(1).phi_deg;
n_elements = numel(patterns);

component_names = sort(fieldnames(patterns(1).components));
element_pattern_stacks = struct();
for k = 1:numel(component_names)
    element_pattern_stacks.(component_names{k}) = stack_component(patterns, component_names{k});
end
ep_reference = element_pattern_stacks.(component_names{1});
directives = config.directives;
polarization = get_field(config, 'polarization', component_names{1});

% ── Weights ───────────────────────────────────────────────────────
if ~isempty(weights_csv_path)
    T = readtable(weights_csv_path);
    amp = T.amplitude;
    phase_deg = T.phase_deg;
    if numel(amp) ~= n_elements
        error('manual_weights_render:CountMismatch', ...
            'CSV has %d rows but %d elements loaded.', numel(amp), n_elements);
    end
    weights_complex = amp(:) .* exp(1i * deg2rad(phase_deg(:)));
else
    weights_complex = ones(n_elements, 1);
end

% ── Power-normalise + array factors (mirrors manual_weights._recompute) ──
w_norm = weights_complex / sqrt(max(sum(abs(weights_complex) .^ 2), 1e-30));

af_components = struct();
for k = 1:numel(component_names)
    af_components.(component_names{k}) = compute_array_factor(w_norm, element_pattern_stacks.(component_names{k}));
end

% Spherical total radiated power (CST partial-directivity normaliser).
theta_rad  = deg2rad(theta_deg(:));
if numel(theta_rad) > 1, dtheta = mean(diff(theta_rad)); else, dtheta = pi; end
if numel(phi_deg) > 1,   dphi = deg2rad(mean(diff(phi_deg(:)))); else, dphi = 2 * pi; end
sin_t   = sin(theta_rad);
p_total = 0;
for k = 1:numel(component_names)
    p_total = p_total + spherical_power(af_components.(component_names{k}), sin_t, dtheta, dphi);
end

if strcmpi(polarization, 'total')
    power_linear_grid = zeros(size(theta_deg, 1), numel(phi_deg));
    for k = 1:numel(component_names)
        power_linear_grid = power_linear_grid + abs(af_components.(component_names{k})) .^ 2;
    end
    metrics_af = sqrt(power_linear_grid);
else
    matches = component_names(strcmpi(component_names, polarization));
    if isempty(matches)
        error('manual_weights_render:BadPolarization', ...
            'Polarization ''%s'' not found. Available components: {%s} (or ''total'').', ...
            polarization, strjoin(component_names, ', '));
    end
    metrics_af = af_components.(matches{1});
end

dbi_grid = compute_directivity_dbi_grid(metrics_af, theta_deg, phi_deg, p_total);
metrics  = evaluate_metrics(ep_reference, theta_deg, phi_deg, weights_complex, ...
    directives, [], metrics_af, p_total);

% ── Heatmap (dBi) with directive overlays ─────────────────────────
fig = figure('Visible', 'off', 'Position', [100 100 900 520]);
ax  = axes(fig); hold(ax, 'on');
pcolor(ax, phi_deg, theta_deg, dbi_grid); shading(ax, 'flat');
colormap(ax, 'jet'); cb = colorbar(ax); cb.Label.String = 'dBi (absolute directivity)';
set(ax, 'YDir', 'reverse'); xlim(ax, [0 360]); ylim(ax, [0 180]);
xlabel(ax, 'Azimuth \phi (deg)'); ylabel(ax, 'Elevation \theta (deg)');
title(ax, sprintf('Array Factor - %s (pol: %s)', basename(patterns_dir), polarization));

phys_masks = build_directive_physical_masks(theta_deg, phi_deg, directives);
for k = 1:numel(directives)
    d = directives{k};
    if strcmp(d.type, 'peak'), color = [0 0.55 0]; else, color = [0.85 0 0]; end
    m = phys_masks{k};
    if any(m(:))
        contour(ax, phi_deg, theta_deg, double(m), [0.5 0.5], 'Color', color, 'LineWidth', 2);
    end
    plot(ax, get_field(d, 'phi', 0.0), d.theta, '+', 'Color', color, 'MarkerSize', 12, 'LineWidth', 2);
end
exportgraphics(fig, output_path, 'Resolution', 150);
close(fig);

% ── Summary ───────────────────────────────────────────────────────
fprintf('Manual render: %d elements, pol=%s\n', n_elements, polarization);
fprintf('  Global peak:  %.2f dBi @ theta=%.1f phi=%.1f\n', ...
    metrics.global_peak_dbi, metrics.global_peak_theta_deg, metrics.global_peak_phi_deg);
fprintf('  3 dB HPBW:    theta %.1f deg, phi %.1f deg\n', metrics.hpbw_theta_deg, metrics.hpbw_phi_deg);
for k = 1:numel(metrics.directive_metrics)
    dm = metrics.directive_metrics(k);
    if strcmp(dm.type, 'peak')
        fprintf('  Peak @ theta=%.1f phi=%.1f: %.2f dBi\n', dm.theta_deg, dm.phi_deg, dm.gain_dbi);
    else
        fprintf('  Null @ theta=%.1f phi=%.1f: %.2f dBi (depth %.1f dB)\n', ...
            dm.theta_deg, dm.phi_deg, dm.gain_dbi, dm.null_depth_db);
    end
end
fprintf('Saved heatmap -> %s\n', output_path);
end


% ────────────────────────── HELPERS ───────────────────────────────

function p = spherical_power(af, sin_t, dtheta, dphi)
p = sum(abs(af) .^ 2 .* sin_t, 'all') * dtheta * dphi;
end


function val = get_field(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end


function b = basename(p)
parts = strsplit(strip_trailing_sep(p), filesep);
b = parts{end};
end


function p = strip_trailing_sep(p)
while ~isempty(p) && (endsWith(p, '/') || endsWith(p, '\'))
    p = p(1:end - 1);
end
end
