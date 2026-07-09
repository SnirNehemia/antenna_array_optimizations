function save_pattern_gif(weights_history, element_patterns_stacked, ...
                          theta_deg, phi_deg, directives, cost_history, ...
                          output_config, output_dir, element_patterns_secondary)
% SAVE_PATTERN_GIF  Animated GIF of the pattern evolving over optimizer iterations.
%
%   Mirrors src/plot/plotter.py:save_pattern_gif. Each frame renders the 2-D
%   heatmap (linear theta) with directive overlays plus a convergence panel.
%   Frames are strided to stay within gif_max_frames. Written with imwrite.
%
%   Inputs:
%       weights_history            : cell array of complex weight vectors per iter.
%       element_patterns_stacked   : (N_elements, N_theta, N_phi) complex.
%       theta_deg, phi_deg         : angle grids [deg].
%       directives                 : cell array of directive structs.
%       cost_history               : cost-per-iteration vector.
%       output_config              : struct (reads gif_max_frames, plot_dynamic_range_db).
%       output_dir                 : directory for the GIF.
%       element_patterns_secondary : optional cross-pol stack ([] -> copol only).
%
%   Part of: Antenna Array Pattern Optimization Tool.

if nargin < 9
    element_patterns_secondary = [];
end
if isempty(weights_history)
    return;
end

GIF_FILENAME = 'pattern_evolution.gif';
MAX_GIF_FRAMES = 100;
GIF_FPS = 10;
DEFAULT_DYNAMIC_RANGE_DB = 40;

theta_deg = theta_deg(:);
phi_deg   = phi_deg(:);

max_frames = cfg_get(output_config, 'gif_max_frames', MAX_GIF_FRAMES);
dyn        = cfg_get(output_config, 'plot_dynamic_range_db', DEFAULT_DYNAMIC_RANGE_DB);

n_iters = numel(weights_history);
stride  = max(1, floor(n_iters / max_frames));
frame_indices = 0:stride:(n_iters - 1);          % 0-based iteration indices
if ~ismember(n_iters - 1, frame_indices)
    frame_indices(end + 1) = n_iters - 1;
end

fig = figure('Visible', 'off', 'Position', [100 100 1200 500]);
ax_map  = subplot(1, 2, 1, 'Parent', fig);
ax_cost = subplot(1, 2, 2, 'Parent', fig);

cost_arr = cost_history(:);
plot(ax_cost, 0:(numel(cost_arr) - 1), cost_arr, 'Color', [0.5 0.5 0.5]);
hold(ax_cost, 'on');
cost_dot = plot(ax_cost, NaN, NaN, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel(ax_cost, 'Iteration'); ylabel(ax_cost, 'Cost J');
title(ax_cost, 'Convergence'); grid(ax_cost, 'on');

phys_masks = build_directive_physical_masks(theta_deg, phi_deg, directives);
outpath = fullfile(output_dir, GIF_FILENAME);
delay   = 1 / GIF_FPS;
first   = true;

for fi = frame_indices
    w  = weights_history{fi + 1};
    af = compute_array_factor(w, element_patterns_stacked);
    if ~isempty(element_patterns_secondary)
        af2   = compute_array_factor(w, element_patterns_secondary);
        power = abs(af) .^ 2 + abs(af2) .^ 2;
    else
        power = abs(af) .^ 2;
    end
    power_db  = 10 * log10(max(power, 1e-30));
    grid_norm = max(power_db - max(power_db, [], 'all'), -dyn);

    cla(ax_map); hold(ax_map, 'on');
    pcolor(ax_map, phi_deg, theta_deg, grid_norm); shading(ax_map, 'flat');
    colormap(ax_map, 'jet'); caxis(ax_map, [-dyn 0]);  % [MATLAB] caxis (clim is R2022a+)
    set(ax_map, 'YDir', 'reverse'); xlim(ax_map, [0 360]); ylim(ax_map, [0 180]);
    xlabel(ax_map, 'Phi (deg)'); ylabel(ax_map, 'Elevation \theta (deg)');

    for k = 1:numel(directives)
        d = directives{k};
        if strcmp(d.type, 'peak'), color = [0 0.55 0]; else, color = [0.85 0 0]; end
        m = phys_masks{k};
        if any(m(:))
            contour(ax_map, phi_deg, theta_deg, double(m), [0.5 0.5], 'Color', color, 'LineWidth', 2);
        end
        plot(ax_map, get_directive_field(d, 'phi', 0.0), d.theta, '+', ...
            'Color', color, 'MarkerSize', 10, 'LineWidth', 2);
    end

    if fi < numel(cost_history), cv = cost_history(fi + 1); else, cv = NaN; end
    title(ax_map, sprintf('Iter %d   J = %.4f', fi, cv));
    set(cost_dot, 'XData', fi, 'YData', cv);
    drawnow;

    frame = getframe(fig);
    [A, map] = rgb2ind(frame2im(frame), 256);
    if first
        imwrite(A, map, outpath, 'gif', 'LoopCount', Inf, 'DelayTime', delay);
        first = false;
    else
        imwrite(A, map, outpath, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
    end
end
close(fig);
end


function val = cfg_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end
