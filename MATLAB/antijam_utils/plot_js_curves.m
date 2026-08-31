function plot_js_curves(sigma_s_db, jn_ratio_db, panels, series_labels, ...
                        fig_title, out_path)
% PLOT_JS_CURVES  Amplitude-plane metrics collapsed onto the J/S axis.
%
%   PLOT_JS_CURVES(sigma_s_db, jn_ratio_db, panels, series_labels, fig_title, out_path)
%
%   The companion to PLOT_AMPLITUDE_HEATMAPS, answering the question the
%   heatmaps raise but cannot settle: how much of the (sigma_s, J/N) plane is
%   really just a function of the jammer-to-SIGNAL ratio J/S = J/N - sigma_s?
%   Every map cell is replotted as a point at x = J/S, colored by its
%   sigma_s_db row. If a metric is governed by J/S alone the rows collapse onto
%   one curve; where they fan out, sigma_s matters independently and the
%   heatmap is the honest picture. Reading a 16x16 grid off a heatmap cannot
%   answer that; a curve plot answers it at a glance.
%
%   Each panel overlays one line per (sigma_s row) x (series). Series are the
%   things being compared — typically {'fixed loading', 'adaptive loading',
%   'oracle'} — and are separated by LINE STYLE, while sigma_s is separated by
%   COLOR (parula, with a colorbar). Two visual channels, two variables, no
%   40-entry legend.
%
%   Inputs:
%       sigma_s_db    : (1 x n_sigma) desired-signal power grid [dB].
%       jn_ratio_db   : (1 x n_jn) jammer-to-noise ratio grid [dB].
%       panels        : struct array, one per panel, with fields
%                       data   : (1 x n_series) cell of (n_sigma x n_jn) maps,
%                                in the same order as series_labels.
%                       title  : char, panel title.
%                       ylabel : char, y-axis label (metric + unit).
%                       yref   : optional scalar; a dashed black reference line
%                                (e.g. the availability floor, or 0 dB for a
%                                loss metric). [] or absent = none.
%       series_labels : (1 x n_series) cell of char, legend names.
%       fig_title     : char, figure super-title.
%       out_path      : PNG destination path.
%
%   R2020a: parula/colorbar only, no clim(), exportgraphics guarded.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P11].

n_sigma   = numel(sigma_s_db);
n_series  = numel(series_labels);
n_panels  = numel(panels);
% Line styles carry the series. More than four series being compared on one
% axis is unreadable regardless of encoding, so that is an error, not a wrap.
styles = {'-', '--', ':', '-.'};
if n_series > numel(styles)
    error('plot_js_curves:TooManySeries', ...
        'plot_js_curves separates series by line style and has %d styles; got %d series.', ...
        numel(styles), n_series);
end

[JJ, SS] = meshgrid(jn_ratio_db, sigma_s_db);
JS = JJ - SS;                              % (n_sigma x n_jn) J/S for every cell

cmap = parula(max(n_sigma, 2));
n_rows = floor(sqrt(n_panels));
n_cols = ceil(n_panels / n_rows);

fig = figure('Visible', 'off', 'Color', 'w', 'Renderer', 'painters', ...
    'Position', [40 40 min(560 * n_cols, 1800), min(430 * n_rows, 1200)]);

for ip = 1:n_panels
    p  = panels(ip);
    ax = subplot(n_rows, n_cols, ip);
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    hold(ax, 'on'); grid(ax, 'on');

    for is = 1:n_series
        m = p.data{is};
        for iy = 1:n_sigma
            x = JS(iy, :);
            y = m(iy, :);
            % Sort by J/S so the line is drawn monotonically; within a row the
            % grid is already ordered, but sorting keeps this correct if a
            % caller ever passes a non-monotonic jn grid.
            [x, ord] = sort(x);
            plot(ax, x, y(ord), styles{is}, 'Color', cmap(iy, :), ...
                'LineWidth', 1.1, 'HandleVisibility', 'off');
        end
    end

    if isfield(p, 'yref') && ~isempty(p.yref)
        xl = [min(JS(:)), max(JS(:))];
        plot(ax, xl, [p.yref p.yref], 'k--', 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end

    xlim(ax, [min(JS(:)), max(JS(:))]);
    xlabel(ax, 'jammer-to-signal ratio  J/S = jn\_ratio\_db - sigma\_s\_db  [dB]');
    ylabel(ax, p.ylabel);
    title(ax, p.title, 'Color', 'k', 'Interpreter', 'none');

    % Style legend: one invisible black proxy line per series, so the legend
    % shows the styles without also showing 16 colors of each.
    if ip == 1 && n_series > 1
        h = gobjects(1, n_series);
        for is = 1:n_series
            h(is) = plot(ax, NaN, NaN, styles{is}, 'Color', 'k', 'LineWidth', 1.2);
        end
        legend(ax, h, series_labels, 'Location', 'best', 'AutoUpdate', 'off', ...
            'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6], ...
            'Interpreter', 'none');
    end

    % Colorbar = the sigma_s axis. caxis is set to the grid's own range so the
    % ticks are real dB values, not colormap indices.
    colormap(ax, cmap);
    if n_sigma > 1
        caxis(ax, [min(sigma_s_db), max(sigma_s_db)]);   % [R2020a] no clim()
        cb = colorbar(ax);
        ylabel(cb, 'desired-signal power  sigma\_s\_db  [dB]', 'Color', 'k');
        set(cb, 'Color', 'k');
    end
end

if exist('sgtitle', 'file')
    sgtitle(fig, fig_title, 'FontWeight', 'bold', 'Color', 'k', ...
        'Interpreter', 'none');
end
exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function exportgraphics_compat(fig, out_path)
% exportgraphics exists from R2020a; fall back to print on older releases —
% and also when exportgraphics itself fails (transient graphics-driver loss).
if exist('exportgraphics', 'file')
    try
        exportgraphics(fig, out_path, 'Resolution', 150);
        return
    catch err
        warning('plot_js_curves:ExportFallback', ...
            'exportgraphics failed for %s (%s); retrying with print.', ...
            out_path, err.message);
    end
end
print(fig, out_path, '-dpng', '-r150');
end
