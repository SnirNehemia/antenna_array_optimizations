function plot_amplitude_heatmaps(sigma_s_db, jn_ratio_db, panels, fig_title, out_path)
% PLOT_AMPLITUDE_HEATMAPS  Metric heatmaps over the (signal, jammer) amplitude plane.
%
%   PLOT_AMPLITUDE_HEATMAPS(sigma_s_db, jn_ratio_db, panels, fig_title, out_path)
%
%   Renders one imagesc panel per performance metric on a shared amplitude
%   plane: desired-signal power on the y axis, jammer-to-noise ratio on the x
%   axis, both in dB relative to the engine's sigma_n^2 = 1 per-element noise
%   floor (sim_engine_init). Because both axes are referenced to the same noise
%   floor, the jammer-to-SIGNAL ratio J/S = jn_ratio_db - sigma_s_db is a family
%   of 45-degree lines; those are overlaid as labelled contours so the two
%   distinct failure regimes (signal-limited at low sigma_s, jammer-limited at
%   high J/S) can be read off directly.
%
%   Inputs:
%       sigma_s_db  : (1 x n_sigma) desired-signal power grid [dB], y axis.
%                     Must be uniformly spaced (imagesc assumption).
%       jn_ratio_db : (1 x n_jn) jammer-to-noise ratio grid [dB], x axis.
%                     Must be uniformly spaced.
%       panels      : struct array, one element per panel, with fields
%                     map        : (n_sigma x n_jn) metric values; NaN cells are
%                                  drawn in the axes background color.
%                     title      : char, panel title.
%                     cbar_label : char, colorbar label (metric + unit).
%                     style      : 'sequential' (parula) | 'diverging' (blue-
%                                  white-red, color scale forced symmetric about
%                                  zero so the sign of a difference map reads
%                                  correctly) | 'categorical' (small integer
%                                  codes 0..n-1, one flat color per code, ticked
%                                  by name on the colorbar via cat_labels).
%                     clim       : optional [lo hi] color limits; [] or absent
%                                  = auto from the map's own finite range.
%                     clim_floor : optional [lo hi] MINIMUM span the auto range
%                                  must cover, i.e. the final range is
%                                  [min(lo_auto, lo), max(hi_auto, hi)]. This is
%                                  what stops a physically-uninteresting metric
%                                  from being auto-stretched into a screaming
%                                  full-scale map: the oracle's dead time is
%                                  0.0-0.1 s out of a 60 s run, but auto-scaling
%                                  to [0, 0.05] painted it saturated yellow and
%                                  read as total failure. Ignored when clim is
%                                  given explicitly.
%                     cat_labels : cell array of names for 'categorical' style,
%                                  one per code 0..numel-1.
%       fig_title   : char, figure super-title.
%       out_path    : PNG destination path.
%
%   Cell values are printed on top of each cell when the grid is small enough
%   (<= 100 cells) to stay legible — a sweep this coarse is read as a table as
%   much as a picture.
%
%   R2020a: uses caxis (not clim), a hand-built diverging colormap (no toolbox),
%   and falls back to print when exportgraphics is unavailable.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6, P9].

n_panels = numel(panels);
n_rows   = floor(sqrt(n_panels));
n_cols   = ceil(n_panels / n_rows);

% 'painters' (software vector renderer) rather than the default opengl: these
% panels are flat 2-D patches and text, so painters loses nothing visually, and
% on the [R2020a] compatibility target it keeps a long batch of figures off the
% hardware path. NOTE this is a no-op on recent releases (the property is
% deprecated there) and it does NOT rescue a session that has already lost
% graphics acceleration mid-run — once that happens even print() fails on a
% trivial plot and only restarting MATLAB helps. exportgraphics_compat below
% degrades as far as it can; the caller writes its data before plotting.
fig = figure('Visible', 'off', 'Color', 'w', 'Renderer', 'painters', ...
    'Position', [40 40 min(520 * n_cols, 1800), min(420 * n_rows, 1200)]);

for i = 1:n_panels
    p  = panels(i);
    ax = subplot(n_rows, n_cols, i);
    draw_panel(ax, sigma_s_db, jn_ratio_db, p);
end

% sgtitle exists from R2018b; guard anyway so the figure still renders without it.
if exist('sgtitle', 'file')
    sgtitle(fig, fig_title, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'none');
end
exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function draw_panel(ax, sigma_s_db, jn_ratio_db, p)
% One heatmap panel: image + colorbar + J/S contours + cell labels.
map = p.map;
p.x_grid = jn_ratio_db;      % stashed for finish_panel (shared by both styles)
p.y_grid = sigma_s_db;
imagesc(ax, jn_ratio_db, sigma_s_db, map);
set(ax, 'YDir', 'normal', 'Color', [0.85 0.85 0.85], ...
    'XColor', 'k', 'YColor', 'k', 'Layer', 'top');
hold(ax, 'on');

finite = map(isfinite(map));
if isfield(p, 'clim') && ~isempty(p.clim)
    lo = p.clim(1); hi = p.clim(2);
elseif isempty(finite)
    lo = 0; hi = 1;                       % all-NaN panel: any valid range
else
    lo = min(finite); hi = max(finite);
    % Widen to the caller's declared minimum span, so a metric that never left
    % its floor is drawn AS a floor rather than auto-stretched to full scale.
    if isfield(p, 'clim_floor') && ~isempty(p.clim_floor)
        lo = min(lo, p.clim_floor(1));
        hi = max(hi, p.clim_floor(2));
    end
end
if strcmp(p.style, 'categorical')
    draw_categorical(ax, p);
    return
end
if strcmp(p.style, 'diverging')
    % Force symmetry about zero so "which side of zero" is the visual message.
    a  = max(abs([lo, hi]));
    if a == 0, a = 1; end
    lo = -a; hi = a;
    colormap(ax, diverging_colormap(256));
else
    colormap(ax, parula(256));
end
if hi <= lo, hi = lo + 1; end             % degenerate (constant) map
caxis(ax, [lo hi]);                       % [R2020a] clim() does not exist yet

cb = colorbar(ax);
ylabel(cb, p.cbar_label, 'Color', 'k');
set(cb, 'Color', 'k');

if numel(map) <= 100
    label_cells(ax, sigma_s_db, jn_ratio_db, map, lo, hi);
end
finish_panel(ax, p);
end


function finish_panel(ax, p)
% J/S contours, grid-aligned ticks, labels — shared by both panel styles.
overlay_js_contours(ax, p.y_grid, p.x_grid);
% Ticks exactly on the sampled grid — every cell is a simulated point, not an
% interpolation, and the axes should not suggest otherwise. Thinned when the
% grid is fine enough that every-cell ticks would collide.
set(ax, 'XTick', thin_ticks(p.x_grid), 'YTick', thin_ticks(p.y_grid), ...
    'TickDir', 'out');
xlabel(ax, 'jammer-to-noise ratio  jn\_ratio\_db  [dB]');
ylabel(ax, 'desired-signal power  sigma\_s\_db  [dB]');
title(ax, p.title, 'Color', 'k', 'Interpreter', 'none');
end


function t = thin_ticks(v)
% Every grid value while they fit; otherwise every k-th, keeping both ends.
max_ticks = 9;
if numel(v) <= max_ticks
    t = v;
    return
end
k = ceil(numel(v) / max_ticks);
t = unique([v(1:k:end), v(end)]);
end


function draw_categorical(ax, p)
% Finish a 'categorical' panel: flat color per integer code, named colorbar
% ticks, no cell-value text (the color IS the value, and the names carry it).
n_cat = numel(p.cat_labels);
caxis(ax, [-0.5, n_cat - 0.5]);           % [R2020a] clim() does not exist yet
colormap(ax, categorical_colormap(n_cat));
cb = colorbar(ax, 'Ticks', 0:(n_cat - 1), 'TickLabels', p.cat_labels);
ylabel(cb, p.cbar_label, 'Color', 'k');
set(cb, 'Color', 'k', 'TickLength', 0);
finish_panel(ax, p);
end


function overlay_js_contours(ax, sigma_s_db, jn_ratio_db)
% Constant jammer-to-SIGNAL ratio lines, J/S = jn_ratio_db - sigma_s_db.
% Drawn as a white halo under a black dashed line so they stay readable over
% both ends of every colormap used here.
[JJ, SS] = meshgrid(jn_ratio_db, sigma_s_db);
js = JJ - SS;
levels = -60:10:60;
levels = levels(levels > min(js(:)) & levels < max(js(:)));
if isempty(levels)
    return
end
contour(ax, jn_ratio_db, sigma_s_db, js, levels, ...
    'LineColor', 'w', 'LineWidth', 1.8, 'HandleVisibility', 'off');
contour(ax, jn_ratio_db, sigma_s_db, js, levels, ...
    'LineColor', 'k', 'LineStyle', '--', 'LineWidth', 0.7, 'HandleVisibility', 'off');

% Labels are placed by hand rather than with clabel: clabel's inline mode puts
% them in the middle of each line, directly on top of the per-cell value text,
% and [R2020a] it returns no text handles when given a contour handle with
% automatic placement, so the labels cannot be restyled after the fact. Each
% line is instead labelled where it LEAVES the plotted area, in the half-cell
% margin that imagesc leaves beyond the outermost cell centers.
dx = axis_step(jn_ratio_db);
dy = axis_step(sigma_s_db);
x0 = min(jn_ratio_db); x1 = max(jn_ratio_db);
y1 = max(sigma_s_db);
for L = levels
    x_at_top = y1 + L;                       % the line is sigma = jn - L
    if x_at_top >= x0 && x_at_top <= x1
        xl = x_at_top;  yl = y1 + 0.34 * dy;      % exits through the top edge
    else
        xl = x1 + 0.34 * dx;  yl = x1 - L;        % exits through the right edge
    end
    text(ax, xl, yl, sprintf('J/S %+d', L), 'FontSize', 6, 'Color', 'k', ...
        'BackgroundColor', 'w', 'Margin', 0.5, 'Clipping', 'on', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
end


function label_cells(ax, sigma_s_db, jn_ratio_db, map, lo, hi)
% Print each cell's value, in white or black depending on cell brightness.
for iy = 1:numel(sigma_s_db)
    for ix = 1:numel(jn_ratio_db)
        v = map(iy, ix);
        if ~isfinite(v)
            continue
        end
        frac = (v - lo) / max(hi - lo, eps);
        if frac < 0.45, col = 'w'; else, col = 'k'; end
        text(ax, jn_ratio_db(ix), sigma_s_db(iy), sprintf('%.1f', v), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 7, 'Color', col);
    end
end
end


function d = axis_step(v)
% Grid spacing of a (possibly single-point) uniform axis vector.
if numel(v) > 1
    d = mean(diff(v));
else
    d = 1;
end
end


function cmap = diverging_colormap(n)
% Blue -> white -> red, built from definition (no toolbox colormaps).
half = floor(n / 2);
rest = n - half;
t = linspace(0, 1, half)';                   % blue  -> white
s = linspace(1, 0, rest)';                   % white -> red
cmap = [t, t, ones(half, 1); ones(rest, 1), s, s];
end


function cmap = categorical_colormap(n)
% Flat, well-separated colors for small-integer code maps. Fixed order so a
% given code keeps its color across figures: red (fail), amber, blue, green
% (pass); beyond four it falls back to a hue sweep.
base = [0.84 0.19 0.15; 0.95 0.68 0.18; 0.20 0.44 0.75; 0.24 0.62 0.32];
if n <= size(base, 1)
    cmap = base(1:n, :);
else
    cmap = hsv2rgb([linspace(0, 0.85, n)', 0.55 * ones(n, 1), 0.80 * ones(n, 1)]);
end
end


function exportgraphics_compat(fig, out_path)
% exportgraphics exists from R2020a; fall back to print on older releases —
% and also when exportgraphics itself fails, which it does on a transient
% graphics-driver loss. print takes a different path and generally survives it;
% losing a figure is not worth losing a sweep that took minutes to compute.
if exist('exportgraphics', 'file')
    try
        exportgraphics(fig, out_path, 'Resolution', 150);
        return
    catch err
        warning('plot_amplitude_heatmaps:ExportFallback', ...
            'exportgraphics failed for %s (%s); retrying with print.', ...
            out_path, err.message);
    end
end
print(fig, out_path, '-dpng', '-r150');
end
