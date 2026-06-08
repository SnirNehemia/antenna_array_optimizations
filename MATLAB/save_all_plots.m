function save_all_plots(array_factor_db_grid, theta_deg, phi_deg, ...
                        weights_complex, directives, ...
                        all_cost_histories, best_run_index, all_run_labels, ...
                        output_config, output_dir)
% SAVE_ALL_PLOTS  Save all enabled result plots (headless) to output_dir.
%
%   Mirrors src/plot/plotter.py:save_all_plots. Each save_* flag in
%   output_config gates one figure. array_factor_db_grid must be pre-computed as
%   10*log10(max(|AF|^2,1e-30)). Figures render with Visible='off' and are
%   written via exportgraphics.
%
%   Inputs:
%       array_factor_db_grid : (N_theta x N_phi) dB.
%       theta_deg, phi_deg   : angle grids [deg].
%       weights_complex      : (N_elements x 1) complex.
%       directives           : cell array of directive structs.
%       all_cost_histories   : cell array of per-restart cost vectors.
%       best_run_index       : 1-based index of the best restart.
%       all_run_labels       : cell array of restart labels.
%       output_config        : struct (config.yaml output section, resolved).
%       output_dir           : directory to write PNGs into.
%
%   Ported from src/plot/plotter.py.
%   Part of: Antenna Array Pattern Optimization Tool.

POLAR_FILE = 'pattern_polar.png';
CART_FILE  = 'pattern_cartesian.png';
PROJ_FILE  = 'pattern_2d.png';
WT_FILE    = 'weights.png';
CH_FILE    = 'cost_history.png';
DEFAULT_DYNAMIC_RANGE_DB = 40;

theta_deg = theta_deg(:);
phi_deg   = phi_deg(:);

output_config    = resolve_output_config(output_config, directives);
dynamic_range_db = cfg_get(output_config, 'plot_dynamic_range_db', DEFAULT_DYNAMIC_RANGE_DB);

if cfg_get(output_config, 'save_polar_plot', false)
    save_polar(array_factor_db_grid, theta_deg, phi_deg, directives, ...
        output_config, dynamic_range_db, fullfile(output_dir, POLAR_FILE));
end
if cfg_get(output_config, 'save_cartesian_plot', false)
    save_cartesian(array_factor_db_grid, theta_deg, phi_deg, directives, ...
        output_config, dynamic_range_db, fullfile(output_dir, CART_FILE));
end
if cfg_get(output_config, 'save_2d_projection_plot', false)
    equal_area = logical(cfg_get(output_config, 'plot_equal_area', true));
    save_projection(array_factor_db_grid, theta_deg, phi_deg, directives, ...
        dynamic_range_db, fullfile(output_dir, PROJ_FILE), equal_area);
end
if cfg_get(output_config, 'save_weight_plots', false)
    save_weights(weights_complex, fullfile(output_dir, WT_FILE));
end
if cfg_get(output_config, 'save_cost_history_plot', false)
    save_cost_history(all_cost_histories, best_run_index, all_run_labels, ...
        fullfile(output_dir, CH_FILE));
end
end


% ────────────────────────── CONFIG HELPERS ────────────────────────

function val = cfg_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end


function oc = resolve_output_config(oc, directives)
if ~isfield(oc, 'plot_cut_type') || ~strcmp(oc.plot_cut_type, 'auto')
    return;
end
target = [];
for k = 1:numel(directives)
    if strcmp(directives{k}.type, 'peak')
        target = directives{k};
        break;
    end
end
if isempty(target) && ~isempty(directives)
    target = directives{1};
end
oc.plot_cut_type = 'theta_cut';
if ~isempty(target)
    oc.plot_phi_deg = get_directive_field(target, 'phi', 0.0);
end
end


function [center_deg, lo_deg, hi_deg, color] = directive_window_bounds(d, oc)
sym_width = get_directive_field(d, 'width', []);
theta_hw  = get_directive_field(d, 'theta_width', sym_width) / 2;
phi_hw    = get_directive_field(d, 'phi_width',   sym_width) / 2;
if strcmp(d.type, 'peak'), color = [0 0.55 0]; else, color = [0.85 0 0]; end
if strcmp(oc.plot_cut_type, 'theta_cut')
    center_deg = d.theta; hw = theta_hw;
else
    center_deg = get_directive_field(d, 'phi', 0.0); hw = phi_hw;
end
lo_deg = center_deg - hw;
hi_deg = center_deg + hw;
end


function [visible, on_front] = directive_on_cut(d, oc)
sym_width = get_directive_field(d, 'width', []);
theta_hw  = get_directive_field(d, 'theta_width', sym_width) / 2;
phi_hw    = get_directive_field(d, 'phi_width',   sym_width) / 2;
d_phi   = get_directive_field(d, 'phi', 0.0);
d_theta = d.theta;
if strcmp(oc.plot_cut_type, 'theta_cut')
    fixed_phi = cfg_get(oc, 'plot_phi_deg', 0.0);
    back_phi  = mod(fixed_phi + 180, 360);
    delta_front = abs(mod(d_phi - fixed_phi + 180, 360) - 180);
    delta_back  = abs(mod(d_phi - back_phi  + 180, 360) - 180);
    on_front = delta_front <= phi_hw;
    visible  = on_front || (delta_back <= phi_hw);
else
    fixed_theta = cfg_get(oc, 'plot_theta_deg', 90.0);
    visible  = abs(d_theta - fixed_theta) <= theta_hw;
    on_front = true;
end
end


% ────────────────────────── PLOT FUNCTIONS ────────────────────────

function save_polar(grid_db, theta, phi, directives, oc, dyn, outpath)
if strcmp(oc.plot_cut_type, 'theta_cut')
    fixed_phi = cfg_get(oc, 'plot_phi_deg', 0.0);
    fi = nearest_index(phi, fixed_phi);
    bi = nearest_index(phi, mod(fixed_phi + 180, 360));
    front_db = grid_db(:, fi); back_db = grid_db(:, bi);
    angle_full = [-flipud(theta); theta];
    power_full = [flipud(back_db); front_db];
    cut_label  = sprintf('Phi = %.1f deg', fixed_phi);
else
    [angle_full, power_full, xlab] = extract_cut(grid_db, theta, phi, oc);
    cut_label = xlab;
end
% Clamp to dynamic range BEFORE plotting. MATLAB reflects any rho outside
% rlim across the origin (adds pi to its angle), so values below -dyn would
% appear as phantom lobes on the opposite side.
power_norm = max(power_full - max(power_full), -dyn);

fig = figure('Visible', 'off', 'Position', [100 100 600 600]);
pax = polaraxes(fig);
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
% rlim must be set BEFORE any polarplot call. If set after, MATLAB reflects
% negative-rho data points across the origin (adds pi to their angle), which
% maps back-hemisphere data onto the front side as a spurious lobe.
rlim(pax, [-dyn 0]);
hold(pax, 'on');
polarplot(pax, deg2rad(angle_full), power_norm, 'LineWidth', 1.5);

for k = 1:numel(directives)
    [vis, on_front] = directive_on_cut(directives{k}, oc);
    if ~vis, continue; end
    [~, lo, hi, color] = directive_window_bounds(directives{k}, oc);
    if strcmp(oc.plot_cut_type, 'theta_cut') && ~on_front
        [lo, hi] = deal(-hi, -lo);   % mirror front window to back half (negate and swap bounds)
    end
    polarplot(pax, deg2rad([lo lo]), [-dyn 0], '--', 'Color', color, 'LineWidth', 1);
    polarplot(pax, deg2rad([hi hi]), [-dyn 0], '--', 'Color', color, 'LineWidth', 1);
end
title(pax, sprintf('Array Pattern - %s', cut_label));
exportgraphics(fig, outpath, 'Resolution', 150);
close(fig);
end


function save_cartesian(grid_db, theta, phi, directives, oc, dyn, outpath)
if strcmp(oc.plot_cut_type, 'theta_cut')
    fixed_phi = cfg_get(oc, 'plot_phi_deg', 0.0);
    fi = nearest_index(phi, fixed_phi);
    bi = nearest_index(phi, mod(fixed_phi + 180, 360));
    front_db = grid_db(:, fi); back_db = grid_db(:, bi);
    angle_full = [-flipud(theta); theta];
    power_full = [flipud(back_db); front_db];
    xlab = sprintf('Theta (deg)  [left: Phi=%.1f | right: Phi=%.1f]', ...
        mod(fixed_phi + 180, 360), fixed_phi);
else
    [angle_full, power_full, xlab] = extract_cut(grid_db, theta, phi, oc);
end
power_norm = power_full - max(power_full);

fig = figure('Visible', 'off', 'Position', [100 100 800 400]);
ax  = axes(fig); hold(ax, 'on');
plot(ax, angle_full, power_norm, 'LineWidth', 1.5);

for k = 1:numel(directives)
    [vis, on_front] = directive_on_cut(directives{k}, oc);
    if ~vis, continue; end
    [center, lo, hi, color] = directive_window_bounds(directives{k}, oc);
    if strcmp(oc.plot_cut_type, 'theta_cut') && ~on_front
        shade_band(ax, -hi, -lo, dyn, color); xline(ax, -center, '--', 'Color', color);
    else
        shade_band(ax, lo, hi, dyn, color);   xline(ax, center, '--', 'Color', color);
    end
end

if strcmp(oc.plot_cut_type, 'theta_cut'), xlim(ax, [-180 180]); end
ylim(ax, [-dyn 0]);
xlabel(ax, xlab); ylabel(ax, 'Power (dB, normalized to peak)');
title(ax, 'Array Pattern - Cartesian'); grid(ax, 'on');
exportgraphics(fig, outpath, 'Resolution', 150);
close(fig);
end


function shade_band(ax, lo, hi, dyn, color)
patch(ax, [lo hi hi lo], [-dyn -dyn 0 0], color, ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');
end


function save_projection(grid_db, theta, phi, directives, dyn, outpath, equal_area)
grid_norm = max(grid_db - max(grid_db, [], 'all'), -dyn);

if equal_area
    y_coords = cos(deg2rad(theta));
    y_label  = 'Elevation \theta [equal-area cos\theta]';
else
    y_coords = theta;
    y_label  = 'Elevation \theta (deg)';
end

fig = figure('Visible', 'off', 'Position', [100 100 900 500]);
ax  = axes(fig); hold(ax, 'on');
% pcolor handles non-uniform y (equal-area) and uniform alike.
pcolor(ax, phi, y_coords, grid_norm);
shading(ax, 'flat');
colormap(ax, 'jet'); clim(ax, [-dyn 0]);
cb = colorbar(ax); cb.Label.String = 'dB (normalized to peak)';

phys_masks = build_directive_physical_masks(theta, phi, directives);
for k = 1:numel(directives)
    d = directives{k};
    if strcmp(d.type, 'peak'), color = [0 0.55 0]; else, color = [0.85 0 0]; end
    m = phys_masks{k};
    if any(m(:))
        contour(ax, phi, y_coords, double(m), [0.5 0.5], 'Color', color, 'LineWidth', 2);
    end
    d_phi = get_directive_field(d, 'phi', 0.0);
    if equal_area, d_y = cos(deg2rad(d.theta)); else, d_y = d.theta; end
    plot(ax, d_phi, d_y, '+', 'Color', color, 'MarkerSize', 12, 'LineWidth', 2);
end

xlabel(ax, 'Phi (deg)'); ylabel(ax, y_label);
title(ax, 'Array Pattern - 2D Projection');
xlim(ax, [0 360]);
if equal_area, ylim(ax, [-1 1]); else, set(ax, 'YDir', 'reverse'); ylim(ax, [0 180]); end
exportgraphics(fig, outpath, 'Resolution', 150);
close(fig);
end


function save_weights(weights_complex, outpath)
n = numel(weights_complex);
idx = (0:n - 1).';
amp = abs(weights_complex(:));
phase_deg = angle(weights_complex(:)) * 180 / pi;

fig = figure('Visible', 'off', 'Position', [100 100 800 500]);
ax1 = subplot(2, 1, 1, 'Parent', fig);
bar(ax1, idx, amp, 0.7);
xlabel(ax1, 'Element index'); ylabel(ax1, 'Amplitude |w_n|');
title(ax1, 'Weight Amplitudes'); grid(ax1, 'on');

ax2 = subplot(2, 1, 2, 'Parent', fig);
bar(ax2, idx, phase_deg, 0.7);
xlabel(ax2, 'Element index'); ylabel(ax2, 'Phase \angle w_n (deg)');
title(ax2, 'Weight Phases'); ylim(ax2, [-180 180]); grid(ax2, 'on');

exportgraphics(fig, outpath, 'Resolution', 150);
close(fig);
end


function save_cost_history(all_cost_histories, best_run_index, all_run_labels, outpath)
fig = figure('Visible', 'off', 'Position', [100 100 800 400]);
ax  = axes(fig); hold(ax, 'on');

for i = 1:numel(all_cost_histories)
    hist  = all_cost_histories{i};
    iters = 0:(numel(hist) - 1);
    if i == best_run_index
        lw = 2; col = [0 0.45 0.74];
    else
        lw = 0.7; col = [0.5 0.5 0.5];
    end
    plot(ax, iters, hist, 'LineWidth', lw, 'Color', col);
end

xlabel(ax, 'Iteration'); ylabel(ax, 'Cost J');
title(ax, 'Optimizer Convergence'); grid(ax, 'on');

if best_run_index >= 1 && best_run_index <= numel(all_run_labels)
    legend(ax, {sprintf('Best - %s', all_run_labels{best_run_index})}, 'Location', 'best');
end
exportgraphics(fig, outpath, 'Resolution', 150);
close(fig);
end


function [ang, pdb, xlab] = extract_cut(grid_db, theta, phi, oc)
if strcmp(oc.plot_cut_type, 'theta_cut')
    fp = cfg_get(oc, 'plot_phi_deg', 0.0);
    pidx = nearest_index(phi, fp);
    pdb  = grid_db(:, pidx);
    ang  = theta(:);
    xlab = sprintf('Theta (deg)  [Phi = %.1f deg]', fp);
else
    ft = cfg_get(oc, 'plot_theta_deg', 90.0);
    tidx = nearest_index(theta, ft);
    pdb  = grid_db(tidx, :).';
    ang  = phi(:);
    xlab = sprintf('Phi (deg)  [Theta = %.1f deg]', ft);
end
end
