function plot_difficulty_scatter(points, fig_title, out_path, pass_pct)
% PLOT_DIFFICULTY_SCATTER  Oracle-tracking score against steering-vector coherence.
%
%   PLOT_DIFFICULTY_SCATTER(points, fig_title, out_path, pass_pct)
%
%   [P12] The blind-spot finder. Every case in a campaign is one marker: x is
%   the measured target/jammer steering coherence (kpi_steering_coherence), y is
%   the headline score. Coherence is the array-independent difficulty
%   coordinate, so a 6-element array and a 4x4 land on the same axis and a
%   single degradation trend should emerge — score high and flat while the
%   jammer is well separated in pattern space, falling as coherence approaches 1.
%
%   The figure is read by looking for what is NOT on that trend. A point sitting
%   well below its neighbours is a case that should have been easy and was not,
%   and its label says which array and which angle to go open a trace for. A
%   coherence SPIKE at a large angular separation (visible in the companion
%   lower panel) is the grating-lobe signature.
%
%   Two panels:
%       top    score vs coherence, one colour per array, one marker per scenario.
%       bottom coherence vs angular separation, same colouring — this is what
%              distinguishes "hard because it is close" from "hard because the
%              array folds this angle onto the target".
%
%   Inputs:
%       points : struct array, one element per case, with fields
%                coh        : steering coherence in [0, 1]. Units: dimensionless.
%                score_pct  : oracle-tracking score. Units: percent.
%                sep_deg    : angular separation target-to-jammer. Units: degrees.
%                array_id   : char, array name (sets colour).
%                scenario_id: char, scenario name (sets marker).
%                label      : char, short annotation drawn next to failing
%                             points only (annotating all of them is unreadable).
%       fig_title : char, figure super-title.
%       out_path  : PNG destination path.
%       pass_pct  : score below which a point is treated as a finding: it is
%                   labelled and ringed. Units: percent.
%
%   R2020a: no toolbox calls; falls back to print when exportgraphics is
%   unavailable.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P12].

% How many failing points get a text label. Rings mark all of them; labels are
% rationed because they overprint.
N_LABELLED = 5;

if isempty(points)
    warning('plot_difficulty_scatter:NoPoints', ...
        'No points to plot; skipping %s.', out_path);
    return
end

arrays    = unique({points.array_id},    'stable');
scenarios = unique({points.scenario_id}, 'stable');
colors  = array_colors(numel(arrays));
markers = {'o', 's', '^', 'd', 'v', 'p', 'h'};

fig = figure('Visible', 'off', 'Color', 'w', 'Renderer', 'painters', ...
    'Position', [40 40 1000 820]);

% ── Top panel: score vs coherence ─────────────────────────────────
ax1 = subplot(2, 1, 1);
hold(ax1, 'on');
legend_handles = gobjects(0);
legend_names   = {};
for ia = 1:numel(arrays)
    for is = 1:numel(scenarios)
        sel = strcmp({points.array_id}, arrays{ia}) & ...
              strcmp({points.scenario_id}, scenarios{is});
        if ~any(sel), continue; end
        mk = markers{mod(is - 1, numel(markers)) + 1};
        h = plot(ax1, [points(sel).coh], [points(sel).score_pct], mk, ...
            'MarkerFaceColor', colors(ia, :), 'MarkerEdgeColor', 'k', ...
            'MarkerSize', 7, 'LineStyle', 'none');
        legend_handles(end + 1) = h;                                     %#ok<AGROW>
        legend_names{end + 1}   = sprintf('%s / %s', arrays{ia}, scenarios{is}); %#ok<AGROW>
    end
end

% Ring every failure, but LABEL only the worst few. The first run of this
% figure labelled all eleven failures and they overprinted each other into an
% unreadable band — the ring already says "this one failed", and the ranked
% list below the axes carries the rest.
fail = [points.score_pct] < pass_pct;
if any(fail)
    plot(ax1, [points(fail).coh], [points(fail).score_pct], 'o', ...
        'MarkerSize', 14, 'MarkerEdgeColor', [0.6 0 0], 'LineWidth', 2, ...
        'LineStyle', 'none');
    idx = find(fail);
    [~, ord] = sort([points(idx).score_pct]);
    idx = idx(ord(1:min(N_LABELLED, numel(ord))));
    for i = 1:numel(idx)
        text(ax1, points(idx(i)).coh, points(idx(i)).score_pct, ...
            ['  ' points(idx(i)).label], 'FontSize', 7, 'Color', [0.4 0 0], ...
            'Interpreter', 'none', 'VerticalAlignment', 'middle', ...
            'Clipping', 'on');
    end
end

yline_compat(ax1, pass_pct, [0.6 0 0], sprintf('pass = %.0f%%', pass_pct));
grid(ax1, 'on');
xlabel(ax1, 'steering coherence |e_s'' e_j| / (|e_s| |e_j|)   [0 = orthogonal, 1 = indistinguishable]');
ylabel(ax1, 'oracle-tracking score [%]');
title(ax1, 'score vs difficulty — points below the trend are the blind spots', ...
    'Color', 'k', 'Interpreter', 'none');
ylim(ax1, [-2 102]);
if numel(legend_handles) <= 12
    lg = legend(ax1, legend_handles, legend_names, 'Location', 'southwest', ...
        'Interpreter', 'none', 'FontSize', 7);
    light_legend(lg);
end
light_axes(ax1);   % after the labels exist, so their Color sticks
hold(ax1, 'off');

% ── Bottom panel: coherence vs angular separation ─────────────────
% Separates "hard because it is close to the target" from "hard because this
% array folds this angle back onto the target" — the latter is a grating lobe
% and shows up here as a coherence spike at a large separation.
ax2 = subplot(2, 1, 2);
hold(ax2, 'on');
for ia = 1:numel(arrays)
    sel = strcmp({points.array_id}, arrays{ia});
    if ~any(sel), continue; end
    [s_sorted, ord] = sort([points(sel).sep_deg]);
    c_all = [points(sel).coh];
    plot(ax2, s_sorted, c_all(ord), '-o', 'Color', colors(ia, :), ...
        'MarkerFaceColor', colors(ia, :), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', arrays{ia});
end
grid(ax2, 'on');
xlabel(ax2, 'angular separation target-to-jammer [deg]');
ylabel(ax2, 'steering coherence');
title(ax2, 'coherence vs separation — a spike at large separation is a grating lobe', ...
    'Color', 'k', 'Interpreter', 'none');
light_legend(legend(ax2, 'Location', 'northeast', 'Interpreter', 'none', ...
    'FontSize', 7));
light_axes(ax2);
hold(ax2, 'off');

if exist('sgtitle', 'file')
    sgtitle(fig, fig_title, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'none');
end
exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function light_legend(lg)
% A legend does not inherit the axes' forced light theme, so it has to be told
% separately or it renders as a black box on a white figure.
set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.5 0.5 0.5]);
end


function light_axes(ax)
% Force the light theme. Recent MATLAB desktops hand a figure dark axes when
% the IDE theme is dark, which produced a black-backgrounded, grey-labelled
% figure on the first run of this suite — unusable in a report and inconsistent
% with every other figure in the milestone, which set these explicitly.
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
    'GridColor', [0.15 0.15 0.15], 'Layer', 'top');
set(get(ax, 'Title'),  'Color', 'k');
set(get(ax, 'XLabel'), 'Color', 'k');
set(get(ax, 'YLabel'), 'Color', 'k');
end


function c = array_colors(n)
% A small qualitative palette, hand-listed so no toolbox colormap is needed on
% the [R2020a] target. Cycles if there are ever more arrays than entries.
base = [0.00 0.45 0.70;      % blue
        0.85 0.37 0.01;      % orange
        0.00 0.62 0.45;      % teal
        0.80 0.47 0.65;      % pink
        0.34 0.34 0.34];     % grey
c = base(mod((0:n - 1), size(base, 1)) + 1, :);
end


function yline_compat(ax, y, color, label)
% [MATLAB] yline exists from R2018b, but drawing the line by hand keeps this
% working unchanged on the R2020a target and gives control over the label
% placement, which yline does not on older releases.
xl = xlim(ax);
plot(ax, xl, [y y], '--', 'Color', color, 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
text(ax, xl(1), y, ['  ' label], 'Color', color, 'FontSize', 8, ...
    'VerticalAlignment', 'bottom', 'Interpreter', 'none');
xlim(ax, xl);
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
        warning('plot_difficulty_scatter:ExportFallback', ...
            'exportgraphics failed for %s (%s); retrying with print.', ...
            out_path, err.message);
    end
end
print(fig, out_path, '-dpng', '-r150');
end
