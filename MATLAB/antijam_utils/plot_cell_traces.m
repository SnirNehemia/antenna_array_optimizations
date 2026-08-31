function plot_cell_traces(scn, aj, traces, dir_ref_dbi, fig_title, out_path)
% PLOT_CELL_TRACES  Full time histories for one cell of the amplitude sweep.
%
%   PLOT_CELL_TRACES(scn, aj, traces, dir_ref_dbi, fig_title, out_path)
%
%   The amplitude sweep collapses each run to six scalars, which is what makes
%   a 256-cell grid readable — and which also hides every transient that
%   produced those scalars. This is the escape hatch: for a handful of
%   deliberately chosen cells (signal-limited, jammer-limited, on the
%   availability cliff, and the high-SNR corner where the beamformer cancels
%   its own signal), re-run the cell and draw what actually happened.
%
%   Three stacked panels on a shared time axis:
%       1. SINR(t) for every algorithm/loading variant, plus the oracle and the
%          sinr_min_db threshold. Time below the threshold is what the sweep's
%          "dead time" cell counts; here you can see whether it is one long
%          outage or a comb of post-toggle dips.
%       2. Directivity toward (theta_s, phi_s) [dBi], against the QUIESCENT
%          beam's directivity as a dashed reference. This is the panel that
%          exposes desired-signal cancellation: a run whose SINR looks healthy
%          while this trace collapses is winning on a numerator it is also
%          destroying, and will fall apart the moment the signal model is wrong.
%       3. Instantaneous oracle gap (oracle SINR - achieved SINR) [dB], which
%          is the sweep's headline metric unrolled in time.
%
%   Jammer-OFF spans are shaded gray in all three panels and scenario events
%   are ticked, so every feature can be attributed to a cause without flipping
%   back to the scenario overview figure.
%
%   Inputs:
%       scn         : struct from sim_scenario (t_s, jammer_on, events).
%       aj          : antijam config section. Required: sinr_min_db.
%       traces      : struct array, one per variant, with fields
%                     label      : char, legend name (e.g. 'lcmv/fixed').
%                     sinr_db    : (1 x T) achieved SINR.
%                     dir_s_dbi  : (1 x T) directivity toward the target.
%                     is_oracle  : logical; the oracle is drawn black/dotted and
%                                  is the reference for panel 3.
%       dir_ref_dbi : scalar, quiescent-beam directivity toward the target
%                     [dBi] — the "no adaptation asked of it" reference for
%                     panel 2. Pass NaN to omit the reference line.
%       fig_title   : char, figure super-title.
%       out_path    : PNG destination path.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P11].

if ~isfield(aj, 'sinr_min_db') || isempty(aj.sinr_min_db)
    error('plot_cell_traces:MissingKey', ...
        'Missing required antijam config key: ''sinr_min_db''.');
end
n_tr = numel(traces);
if n_tr == 0
    error('plot_cell_traces:NoTraces', 'traces is empty — nothing to draw.');
end
i_oracle = find([traces.is_oracle], 1);

t     = scn.t_s;
t_max = t(end);
% Colors are assigned to the NON-oracle traces only (the oracle is always black
% and dotted), from a hand-picked palette rather than lines(): lines() spends a
% slot on the oracle and its third entry is yellow, which is unreadable as a
% thin line on a white background.
cols = trace_colors(n_tr, [traces.is_oracle]);

fig = figure('Visible', 'off', 'Color', 'w', 'Renderer', 'painters', ...
    'Position', [40 40 1100 780]);
ax = gobjects(1, 3);
for i = 1:3
    ax(i) = subplot(3, 1, i);
    set(ax(i), 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    hold(ax(i), 'on'); grid(ax(i), 'on');
end

% ── 1. SINR ───────────────────────────────────────────────────────
all_sinr = cell2mat({traces.sinr_db}');
span = pad_span([min(all_sinr(:)), max(all_sinr(:)), aj.sinr_min_db]);
ylim(ax(1), span);
shade_off_spans(ax(1), t, ~scn.jammer_on, span);
mark_events(ax(1), scn, span);
plot(ax(1), [0 t_max], aj.sinr_min_db * [1 1], 'r--', 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
h = gobjects(1, n_tr);
for i = 1:n_tr
    h(i) = plot_variant(ax(1), t, traces(i).sinr_db, traces(i), cols(i, :));
end
ylabel(ax(1), 'SINR [dB]');
title(ax(1), sprintf(['achieved SINR (red: %.1f dB threshold; gray: jammer ' ...
    'OFF)'], aj.sinr_min_db), 'Color', 'k');
legend(ax(1), h, {traces.label}, 'Location', 'best', 'AutoUpdate', 'off', ...
    'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6], ...
    'Interpreter', 'none');

% ── 2. Directivity toward the target ──────────────────────────────
all_dir = cell2mat({traces.dir_s_dbi}');
ref_vals = all_dir(:);
if isfinite(dir_ref_dbi), ref_vals = [ref_vals; dir_ref_dbi]; end
span = pad_span([min(ref_vals), max(ref_vals)]);
ylim(ax(2), span);
shade_off_spans(ax(2), t, ~scn.jammer_on, span);
mark_events(ax(2), scn, span);
if isfinite(dir_ref_dbi)
    plot(ax(2), [0 t_max], dir_ref_dbi * [1 1], 'k--', 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end
for i = 1:n_tr
    plot_variant(ax(2), t, traces(i).dir_s_dbi, traces(i), cols(i, :));
end
ylabel(ax(2), 'directivity @ target [dBi]');
title(ax(2), sprintf(['directivity toward (\\theta_s, \\phi_s) — dashed: ' ...
    'quiescent beam %.2f dBi. A trace far BELOW it is cancelling the ' ...
    'desired signal.'], dir_ref_dbi), 'Color', 'k');

% ── 3. Instantaneous oracle gap ───────────────────────────────────
if isempty(i_oracle)
    axis(ax(3), 'off');
    text(ax(3), 0.5, 0.5, 'no oracle trace supplied', 'Color', 'k', ...
        'HorizontalAlignment', 'center');
else
    gaps = [];
    for i = 1:n_tr
        if i == i_oracle, continue, end
        gaps = [gaps; traces(i_oracle).sinr_db - traces(i).sinr_db]; %#ok<AGROW>
    end
    span = pad_span([min(gaps(:)), max(gaps(:)), 0]);
    ylim(ax(3), span);
    shade_off_spans(ax(3), t, ~scn.jammer_on, span);
    mark_events(ax(3), scn, span);
    plot(ax(3), [0 t_max], [0 0], 'k--', 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    for i = 1:n_tr
        if i == i_oracle, continue, end
        plot(ax(3), t, traces(i_oracle).sinr_db - traces(i).sinr_db, '-', ...
            'Color', cols(i, :), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    ylabel(ax(3), 'oracle gap [dB]');
    title(ax(3), 'instantaneous oracle gap (oracle SINR - achieved; lower is better)', ...
        'Color', 'k');
end
xlabel(ax(3), 't [s]');

for i = 1:3
    xlim(ax(i), [0 t_max]);
end

if exist('sgtitle', 'file')
    sgtitle(fig, fig_title, 'FontWeight', 'bold', 'Color', 'k', ...
        'Interpreter', 'none');
end
exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function cols = trace_colors(n_tr, is_oracle)
% One row per trace; oracle rows are placeholders (never used — plot_variant
% draws the oracle black). Non-oracle traces take the palette in order.
palette = [0.85 0.33 0.10;    % orange
           0.00 0.45 0.74;    % blue
           0.47 0.67 0.19;    % green
           0.49 0.18 0.56;    % purple
           0.64 0.08 0.18];   % dark red
cols = zeros(n_tr, 3);
k = 0;
for i = 1:n_tr
    if is_oracle(i)
        continue
    end
    k = k + 1;
    cols(i, :) = palette(mod(k - 1, size(palette, 1)) + 1, :);
end
end


function h = plot_variant(ax, t, y, tr, col)
% The oracle is always black and dotted; everything else takes its series color.
if tr.is_oracle
    h = plot(ax, t, y, 'k:', 'LineWidth', 1.4);
else
    h = plot(ax, t, y, '-', 'Color', col, 'LineWidth', 1.2);
end
end


function s = pad_span(vals)
% A finite y span with 8% headroom, robust to a constant trace.
vals = vals(isfinite(vals));
if isempty(vals)
    s = [0 1];
    return
end
lo = min(vals); hi = max(vals);
if hi <= lo, hi = lo + 1; end
pad = 0.08 * (hi - lo);
s = [lo - pad, hi + pad];
end


function shade_off_spans(ax, t_s, off_mask, y_range)
% Gray vertical spans over every contiguous run of `off_mask`. Explicit y
% extents, not +-inf: a patch with huge coordinates drags the axis autoscale.
if ~any(off_mask)
    return
end
d      = diff([false, off_mask, false]);
starts = find(d == 1);
stops  = find(d == -1) - 1;
for i = 1:numel(starts)
    x0 = t_s(starts(i));
    x1 = t_s(stops(i));
    patch(ax, [x0 x1 x1 x0], y_range([1 1 2 2]), [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
end
end


function mark_events(ax, scn, y_range)
% Faint vertical tick at every scenario event (turn_on / power_step / jump...).
for k = 1:numel(scn.events)
    plot(ax, scn.events{k}.t_s * [1 1], y_range, ':', ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
end
end


function exportgraphics_compat(fig, out_path)
% exportgraphics exists from R2020a; fall back to print on older releases —
% and also when exportgraphics itself fails (transient graphics-driver loss).
if exist('exportgraphics', 'file')
    try
        exportgraphics(fig, out_path, 'Resolution', 150);
        return
    catch err
        warning('plot_cell_traces:ExportFallback', ...
            'exportgraphics failed for %s (%s); retrying with print.', ...
            out_path, err.message);
    end
end
print(fig, out_path, '-dpng', '-r150');
end
