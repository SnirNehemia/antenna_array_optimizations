function plot_scenario_overview(scenarios, aj, fig_title, out_path)
% PLOT_SCENARIO_OVERVIEW  What the jammer actually does, before any algorithm runs.
%
%   PLOT_SCENARIO_OVERVIEW(scenarios, aj, fig_title, out_path)
%
%   Draws the GROUND TRUTH timeline of one or more scenarios: where the jammer
%   is, when it transmits, how hard, and how far off the main beam it sits. It
%   consumes only the struct returned by sim_scenario — no closed-loop run, no
%   beamformer, no seeds — so it is the reference picture a results folder needs
%   in order to be readable on its own. Without it the metric heatmaps carry no
%   statement of what they were measured against, and an "ONOFF" label is a
%   promise rather than a plot.
%
%   Four stacked panels, shared time axis, one colored line per scenario:
%       1. theta_j(t)   elevation [deg], with the target elevation as a dashed
%                       reference and the guard cap as a shaded band.
%       2. phi_j(t)     azimuth [deg], same references. Drawn with NaN breaks
%                       inserted at the 0/360 wrap so a wrapping trajectory does
%                       not draw a vertical line across the whole panel.
%       3. jn_ratio_db(t) with the OFF phases masked out (NaN) and shaded, so
%                       duty cycle, toggle period and power steps are all read
%                       off the same axis. This is the panel that answers "is
%                       the jammer on right now".
%       4. angular separation from (theta_s, phi_s) [deg] — the single number
%                       that says how hard the geometry is, with the guard_deg
%                       floor marked. A scenario whose separation approaches
%                       guard_deg is asking the nuller to work inside the main
%                       lobe, which is out of scope for this milestone.
%
%   Event markers (sim_scenario's events cell: 'jump' / 'power_step' /
%   'turn_on' / 'turn_off' / 'freq_hop') are drawn as vertical ticks on the
%   power panel and, for a single-scenario figure, labelled by type. With many
%   scenarios overlaid the labels are dropped and only the ticks remain.
%
%   NOTE the ON/OFF shading and every value here are simulator ground truth,
%   the same arrays fed to sim_engine_step and the oracle. No adapt_/agent_
%   algorithm sees any of it (plan constraint, sim_scenario.m header).
%
%   Inputs:
%       scenarios : struct array OR cell array of structs from sim_scenario.
%                   All must share the same t_s length only if you want them to
%                   line up visually; different durations are fine (the x axis
%                   spans the longest).
%       aj        : antijam config section. Required: theta_s_deg, phi_s_deg,
%                   guard_deg.
%       fig_title : char, figure super-title.
%       out_path  : PNG destination path.
%
%   R2020a: no clim/exportgraphics assumptions beyond the shared compat helper.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P11].

if isstruct(scenarios)
    scenarios = num2cell(scenarios);
end
if ~iscell(scenarios) || isempty(scenarios)
    error('plot_scenario_overview:BadScenarios', ...
        'scenarios must be a non-empty struct array or cell array of sim_scenario structs.');
end
for key = {'theta_s_deg', 'phi_s_deg', 'guard_deg'}
    if ~isfield(aj, key{1}) || isempty(aj.(key{1}))
        error('plot_scenario_overview:MissingKey', ...
            'Missing required antijam config key: ''%s''.', key{1});
    end
end

n_scn  = numel(scenarios);
cols   = lines(max(n_scn, 3));
% Widths decrease with scenario index so that identical trajectories (STATIC,
% ONOFF and FASTONOFF share one jammer position by design) stay individually
% visible as concentric bands, instead of the last one drawn hiding the rest.
lw     = max(2.6 - 0.45 * (0:n_scn - 1), 0.9);
t_max  = max(cellfun(@(s) s.t_s(end), scenarios));
single = (n_scn == 1);

fig = figure('Visible', 'off', 'Color', 'w', 'Renderer', 'painters', ...
    'Position', [40 40 1100 900]);

ax = gobjects(1, 4);
for i = 1:4
    ax(i) = subplot(4, 1, i);
    set(ax(i), 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    hold(ax(i), 'on'); grid(ax(i), 'on');
end

% ── 1. Elevation ──────────────────────────────────────────────────
% Guard band drawn first so the trajectories sit on top of it.
shade_band(ax(1), [0 t_max], aj.theta_s_deg + [-aj.guard_deg, aj.guard_deg]);
yline_compat(ax(1), aj.theta_s_deg, 'k--', [0 t_max]);
for i = 1:n_scn
    plot(ax(1), scenarios{i}.t_s, scenarios{i}.theta_j_deg, '-', ...
        'Color', cols(i, :), 'LineWidth', lw(i));
end
ylabel(ax(1), '\theta_j [deg]');
title(ax(1), sprintf(['jammer elevation (dashed: target \\theta_s = %.0f' ...
    '\\circ; shaded: \\pm%.0f\\circ — but the guard cap is a 2-D cone, so ' ...
    'only the separation panel decides)'], aj.theta_s_deg, aj.guard_deg), ...
    'Color', 'k');

% ── 2. Azimuth ────────────────────────────────────────────────────
shade_band(ax(2), [0 t_max], aj.phi_s_deg + [-aj.guard_deg, aj.guard_deg]);
yline_compat(ax(2), aj.phi_s_deg, 'k--', [0 t_max]);
for i = 1:n_scn
    [tw, pw] = break_at_wrap(scenarios{i}.t_s, scenarios{i}.phi_j_deg);
    plot(ax(2), tw, pw, '-', 'Color', cols(i, :), 'LineWidth', lw(i));
end
ylabel(ax(2), '\phi_j [deg]');
ylim(ax(2), [0 360]);
set(ax(2), 'YTick', 0:90:360);
title(ax(2), sprintf('jammer azimuth (dashed: target \\phi_s = %.0f\\circ)', ...
    aj.phi_s_deg), 'Color', 'k');

% ── 3. Power / on-off ─────────────────────────────────────────────
% OFF phases: shaded span + a NaN break in the power line. Two encodings of the
% same fact on purpose — the shading survives being read at thumbnail size, the
% gap survives being read in grayscale.
% Fix the y span first: the OFF shading and the event ticks have to be drawn to
% finite extents, or the patches drag the autoscale off to their own corners.
all_jn = cell2mat(cellfun(@(s) s.jn_ratio_db(:).', scenarios, 'UniformOutput', false));
jn_lo  = min(all_jn); jn_hi = max(all_jn);
if jn_hi <= jn_lo, jn_hi = jn_lo + 1; end
pad    = 0.15 * (jn_hi - jn_lo);
jn_span = [jn_lo - pad, jn_hi + pad];
ylim(ax(3), jn_span);
for i = 1:n_scn
    s = scenarios{i};
    shade_off_spans(ax(3), s.t_s, ~s.jammer_on, jn_span, 0.6 / n_scn);
    jn = s.jn_ratio_db;
    jn(~s.jammer_on) = NaN;
    plot(ax(3), s.t_s, jn, '-', 'Color', cols(i, :), 'LineWidth', lw(i));
end
mark_events(ax(3), scenarios, cols, single, jn_span);
ylabel(ax(3), 'J/N [dB]');
title(ax(3), 'jammer power (gaps + gray shading = transmitter OFF; ticks = events)', ...
    'Color', 'k');

% ── 4. Angular separation from the target ─────────────────────────
for i = 1:n_scn
    s   = scenarios{i};
    sep = angular_separation_deg(s.theta_j_deg, s.phi_j_deg, ...
        aj.theta_s_deg, aj.phi_s_deg);
    plot(ax(4), s.t_s, sep, '-', 'Color', cols(i, :), 'LineWidth', lw(i));
end
yline_compat(ax(4), aj.guard_deg, 'r--', [0 t_max]);
ylabel(ax(4), 'separation [deg]');
xlabel(ax(4), 't [s]');
title(ax(4), sprintf(['angular separation from the target (red: %.0f\\circ ' ...
    'guard floor — below it the jammer is inside the main beam, out of scope)'], ...
    aj.guard_deg), 'Color', 'k');

for i = 1:4
    xlim(ax(i), [0 t_max]);
end

% One legend for the whole figure, on the top panel (the scenarios are the
% same set of lines in every panel).
labels = cellfun(@(s) s.id, scenarios, 'UniformOutput', false);
if ~single
    legend(ax(1), labels, 'Location', 'best', 'AutoUpdate', 'off', ...
        'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);
end

if exist('sgtitle', 'file')
    sgtitle(fig, fig_title, 'FontWeight', 'bold', 'Color', 'k', ...
        'Interpreter', 'none');
end
exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function shade_band(ax, x_range, y_range)
% A light horizontal band across the full x span (guard cap).
patch(ax, [x_range(1) x_range(2) x_range(2) x_range(1)], ...
    [y_range(1) y_range(1) y_range(2) y_range(2)], [1 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end


function shade_off_spans(ax, t_s, off_mask, y_range, alpha)
% Gray vertical spans over every contiguous run of `off_mask`.
if ~any(off_mask)
    return
end
d     = diff([false, off_mask, false]);
starts = find(d == 1);
stops  = find(d == -1) - 1;
yl = y_range;
for i = 1:numel(starts)
    x0 = t_s(starts(i));
    x1 = t_s(stops(i));
    patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceAlpha', alpha, 'HandleVisibility', 'off');
end
end


function mark_events(ax, scenarios, cols, single, y_range)
% Vertical ticks at every scenario event; typed labels only when there is one
% scenario on the figure AND few enough events to read (a 5 s toggle over 100 s
% produces 19 of them, and 19 rotated labels are a smear, not information).
MAX_LABELS = 8;
for i = 1:numel(scenarios)
    ev = scenarios{i}.events;
    label_them = single && numel(ev) <= MAX_LABELS;
    for k = 1:numel(ev)
        xline_compat(ax, ev{k}.t_s, cols(i, :), y_range);
        if label_them
            text(ax, ev{k}.t_s, y_range(2) - 0.04 * diff(y_range), ...
                [' ' ev{k}.type], 'Rotation', 90, ...
                'FontSize', 6, 'Color', cols(i, :), 'Interpreter', 'none', ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end
    end
end
end


function [t_out, p_out] = break_at_wrap(t_s, phi_deg)
% Insert NaN wherever phi jumps by more than 180 deg (a 0/360 wrap), so the
% plotted line lifts the pen instead of drawing a spurious vertical segment.
jump = [false, abs(diff(phi_deg)) > 180];
t_out = t_s;
p_out = phi_deg;
p_out(jump) = NaN;
end


function yline_compat(ax, y, spec, x_range)
% yline exists from R2018b, but is avoided here for a more basic reason: these
% rules are drawn BEFORE xlim is set, so the span must be passed in explicitly
% rather than read back off the axes.
plot(ax, x_range, [y y], spec, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end


function xline_compat(ax, x, color, y_range)
% Vertical rule over an explicit y span. Explicit, not +-inf, because a patch
% or line with huge coordinates drags the axis autoscale with it.
plot(ax, [x x], y_range, ':', 'Color', color, 'LineWidth', 0.9, ...
    'HandleVisibility', 'off');
end


function exportgraphics_compat(fig, out_path)
% exportgraphics exists from R2020a; fall back to print on older releases —
% and also when exportgraphics itself fails (transient graphics-driver loss).
if exist('exportgraphics', 'file')
    try
        exportgraphics(fig, out_path, 'Resolution', 150);
        return
    catch err
        warning('plot_scenario_overview:ExportFallback', ...
            'exportgraphics failed for %s (%s); retrying with print.', ...
            out_path, err.message);
    end
end
print(fig, out_path, '-dpng', '-r150');
end
