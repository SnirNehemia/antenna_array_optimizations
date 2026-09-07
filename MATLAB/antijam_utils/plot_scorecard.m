function plot_scorecard(row_labels, col_labels, panels, fig_title, out_path, pass_pct)
% PLOT_SCORECARD  Pass/fail scorecard over (array x jammer position), one panel per scenario.
%
%   PLOT_SCORECARD(row_labels, col_labels, panels, fig_title, out_path, pass_pct)
%
%   [P12] The primary "how is the algorithm doing" figure. Deliberately small,
%   dense and readable at a glance: one cell per (array, jammer position), the
%   headline score printed inside it, and a heavy border around anything that
%   failed. It exists to replace the question "which of the nine heatmap panels
%   do I read" with a single picture that has one number in each cell.
%
%   Colour runs red -> amber -> green over a FIXED [0, 100] scale, not the
%   project's usual parula. Two reasons: a scorecard is read as pass/fail rather
%   than as a field, and a fixed scale is what lets two scorecards (two
%   algorithms, or before/after a change) be compared by flipping between them.
%   Auto-scaling would make a uniformly-excellent grid look identical to a
%   uniformly-terrible one.
%
%   Inputs:
%       row_labels : (1 x n_rows) cell of char, y-axis names (arrays).
%       col_labels : (1 x n_cols) cell of char, x-axis names (jammer positions).
%       panels     : struct array, one element per scenario, with fields
%                    map   : (n_rows x n_cols) score values [%]; NaN cells are
%                            drawn in the axes background colour (a case that
%                            failed to run is visibly absent, not silently 0).
%                    title : char, panel title (the scenario id).
%                    note  : optional (n_rows x n_cols) cell of char, a short
%                            second line printed under the score (used for the
%                            steering coherence of that cell). [] or absent to
%                            omit.
%       fig_title  : char, figure super-title.
%       out_path   : PNG destination path.
%       pass_pct   : score at or above which a cell passes; cells below get the
%                    heavy border. Units: percent.
%
%   R2020a: uses caxis (not clim), a hand-built colormap (no toolbox), and
%   falls back to print when exportgraphics is unavailable.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P12].

n_panels = numel(panels);
n_rows_f = floor(sqrt(n_panels));
n_cols_f = ceil(n_panels / n_rows_f);

% 'painters' for the same reason as plot_amplitude_heatmaps: flat 2-D patches
% and text lose nothing on the software renderer, and it keeps a long batch of
% figures off the hardware path on the [R2020a] compatibility target.
% The width floor is not cosmetic: a single-panel scorecard sized at 560 px
% clipped its own super-title on the first run of this suite, and the title is
% where the seed count and the amplitudes are recorded.
fig = figure('Visible', 'off', 'Color', 'w', 'Renderer', 'painters', ...
    'Position', [40 40 max(min(560 * n_cols_f, 1800), 900), ...
                 min(420 * n_rows_f, 1200)]);

for i = 1:n_panels
    ax = subplot(n_rows_f, n_cols_f, i);
    draw_scorecard_panel(ax, row_labels, col_labels, panels(i), pass_pct);
end

if exist('sgtitle', 'file')
    % Wrapped rather than widened. These titles carry the seed count and the
    % amplitudes, and a one-panel scorecard is narrow enough that the string
    % ran off both edges — widening the figure to fit a title just stretches
    % the cells. sgtitle takes a cellstr as multiple lines.
    sgtitle(fig, wrap_title(fig_title, 72), 'FontWeight', 'bold', ...
        'Color', 'k', 'Interpreter', 'none');
end
exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function draw_scorecard_panel(ax, row_labels, col_labels, p, pass_pct)
% One scorecard panel: image + per-cell score text + fail borders.
map = p.map;
[n_rows, n_cols] = size(map);

imagesc(ax, 1:n_cols, 1:n_rows, map);
set(ax, 'YDir', 'reverse', 'Color', [0.85 0.85 0.85], ...
    'XColor', 'k', 'YColor', 'k', 'Layer', 'top', ...
    'XTick', 1:n_cols, 'XTickLabel', col_labels, ...
    'YTick', 1:n_rows, 'YTickLabel', row_labels, ...
    'TickLabelInterpreter', 'none');
caxis(ax, [0 100]);                          % [MATLAB] R2020a: caxis, not clim
colormap(ax, red_amber_green());
hold(ax, 'on');

cb = colorbar(ax);
cb.Label.String = 'oracle-tracking score [%]';
cb.Color = 'k';

has_note = isfield(p, 'note') && ~isempty(p.note);
for ir = 1:n_rows
    for ic = 1:n_cols
        v = map(ir, ic);
        if ~isfinite(v)
            text(ax, ic, ir, 'n/a', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'Color', [0.3 0.3 0.3], ...
                'FontSize', 9);
            continue
        end
        % Text colour picked against the cell, not fixed: the amber midrange is
        % light enough that white text on it is unreadable.
        txt_color = 'k';
        if v < 35, txt_color = 'w'; end
        label = sprintf('%.0f', v);
        if has_note && ~isempty(p.note{ir, ic})
            label = sprintf('%.0f\n%s', v, p.note{ir, ic});
        end
        text(ax, ic, ir, label, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', txt_color, ...
            'FontSize', 9, 'FontWeight', 'bold');
        if v < pass_pct
            % A heavy border rather than a colour change: the colour already
            % encodes the score, and a failing cell must survive being printed
            % in greyscale or looked at by someone who reads red-green poorly.
            rectangle(ax, 'Position', [ic - 0.5, ir - 0.5, 1, 1], ...
                'EdgeColor', [0.6 0 0], 'LineWidth', 3);
        end
    end
end

n_fail = sum(map(:) < pass_pct & isfinite(map(:)));
title(ax, sprintf('%s  (%d/%d below %.0f%%)', p.title, n_fail, ...
    sum(isfinite(map(:))), pass_pct), 'Color', 'k', 'Interpreter', 'none');
xlabel(ax, 'jammer position', 'Color', 'k');
xtickangle(ax, 20);
% Re-assert the light theme after the labels exist: a dark MATLAB desktop hands
% out dark axes, and the colours set at creation do not survive every label
% call. Every other figure in this milestone forces light explicitly.
set(ax, 'Color', [0.85 0.85 0.85], 'XColor', 'k', 'YColor', 'k');
set(get(ax, 'Title'),  'Color', 'k');
set(get(ax, 'XLabel'), 'Color', 'k');
hold(ax, 'off');
end


function lines = wrap_title(txt, width)
% Greedy word wrap to a cellstr, one entry per line. Base MATLAB only — the
% Text Analytics Toolbox is not available on the [R2020a] target.
words = strsplit(txt, ' ');
lines = {};
cur   = '';
for i = 1:numel(words)
    if isempty(cur)
        cand = words{i};
    else
        cand = [cur ' ' words{i}];
    end
    if numel(cand) > width && ~isempty(cur)
        lines{end + 1} = cur;                                         %#ok<AGROW>
        cur = words{i};
    else
        cur = cand;
    end
end
if ~isempty(cur), lines{end + 1} = cur; end
end


function cmap = red_amber_green()
% Hand-built red -> amber -> green ramp. Built by interpolation rather than
% taken from a toolbox because the project targets base [R2020a] with only the
% Optimization Toolbox available.
anchors = [0.70 0.10 0.10;      % 0%   deep red
           0.90 0.45 0.15;      % 33%  orange
           0.95 0.85 0.30;      % 66%  amber
           0.20 0.65 0.30];     % 100% green
n = 256;
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = [interp1(x, anchors(:, 1), xi).', ...
        interp1(x, anchors(:, 2), xi).', ...
        interp1(x, anchors(:, 3), xi).'];
end


function exportgraphics_compat(fig, out_path)
% exportgraphics exists from R2020a; fall back to print on older releases —
% and also when exportgraphics itself fails, which it does on a transient
% graphics-driver loss.
if exist('exportgraphics', 'file')
    try
        exportgraphics(fig, out_path, 'Resolution', 150);
        return
    catch err
        warning('plot_scorecard:ExportFallback', ...
            'exportgraphics failed for %s (%s); retrying with print.', ...
            out_path, err.message);
    end
end
print(fig, out_path, '-dpng', '-r150');
end
