function stats = plot_mode_c_comparison(runs, aj, out_path)
% PLOT_MODE_C_COMPARISON  Side-by-side comparison of Mode C algorithms.
%
%   stats = PLOT_MODE_C_COMPARISON(runs, aj, out_path)
%
%   One figure comparing every algorithm in `runs` on a single scenario:
%     - top panel : achieved SINR vs time for each algorithm, the oracle
%       reference (black dashed), the sinr_min_db threshold (red dashed) and
%       jammer-on shading — the notches at each jammer turn-on are where the
%       methods differ;
%     - three bar panels: SINR availability [%], total 'dead time' [s]
%       (time below threshold), and mean oracle gap [dB] per algorithm.
%
%   Also returns the per-algorithm statistics it plots, so the demo script can
%   print them as a table.
%
%   Inputs:
%       runs     : cell array of run records, each a struct with fields
%                  algorithm (char), run_log (with sinr_db, oracle_sinr_db),
%                  kpi (from kpi_evaluate), scenario (from sim_scenario). The
%                  oracle run, if present, is drawn as the reference line and
%                  excluded from the oracle-gap bars.
%       aj       : antijam config (needs sinr_min_db).
%       out_path : PNG destination path.
%
%   Outputs:
%       stats : struct array (one per run) with fields algorithm, availability
%               (%), dead_time_s, oracle_gap_db, recovery_mean_steps.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

thr = aj.sinr_min_db;
n   = numel(runs);
scn = runs{1}.scenario;
t   = scn.t_s;
dt  = t(2) - t(1);

% ── Per-algorithm statistics ───────────────────────────────────────
stats = struct('algorithm', {}, 'availability', {}, 'dead_time_s', {}, ...
    'oracle_gap_db', {}, 'recovery_mean_steps', {});
labels   = cell(1, n);
avail    = zeros(1, n);
deadtime = zeros(1, n);
gap      = NaN(1, n);
recov    = NaN(1, n);
is_oracle = false(1, n);
for i = 1:n
    r = runs{i};
    sinr = r.run_log.sinr_db;
    labels{i}   = r.algorithm;
    avail(i)    = 100 * mean(sinr >= thr);
    deadtime(i) = dt * sum(sinr < thr);
    is_oracle(i) = strcmp(r.algorithm, 'oracle');
    if ~is_oracle(i)
        gap(i) = r.kpi.oracle_gap_mean_db;
    end
    recov(i) = mean(r.kpi.recovery_steps, 'omitnan');
    stats(i) = struct('algorithm', r.algorithm, 'availability', avail(i), ...
        'dead_time_s', deadtime(i), 'oracle_gap_db', gap(i), ...
        'recovery_mean_steps', recov(i));
end

% Consistent colors: oracle black, others from lines().
cols = lines(max(1, sum(~is_oracle)));
col_of = containers.Map('KeyType', 'char', 'ValueType', 'any');
ci = 0;
for i = 1:n
    if is_oracle(i)
        col_of(labels{i}) = [0 0 0];
    else
        ci = ci + 1;
        col_of(labels{i}) = cols(ci, :);
    end
end

% ── Figure ─────────────────────────────────────────────────────────
fig = figure('Visible', 'off', 'Position', [50 50 1000 700], 'Color', 'w');

% Top: SINR timelines.
axT = subplot(2, 1, 1); hold(axT, 'on'); grid(axT, 'on');
set(axT, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
shade_jammer_on(axT, scn);                         % gray bands while jammer ON
all_sinr = [];
for i = 1:n
    ls = '-'; lw = 1.4;
    if is_oracle(i), ls = '--'; lw = 1.1; end
    plot(axT, t, runs{i}.run_log.sinr_db, ls, 'Color', col_of(labels{i}), ...
        'LineWidth', lw);
    all_sinr = [all_sinr, runs{i}.run_log.sinr_db]; %#ok<AGROW>
end
yline(axT, thr, 'r--', 'LineWidth', 1.2);
% Fix the y-range to the SINR data (the shading patches are intentionally tall).
ylim(axT, [min([all_sinr, thr]) - 2, max(all_sinr) + 3]);
xlim(axT, [t(1), t(end)]);
xlabel(axT, 't [s]'); ylabel(axT, 'SINR [dB]');
title(axT, sprintf('Mode C comparison — scenario %s (jammer on/off)', ...
    strrep(scn.id, '_', '\_')), 'Color', 'k');
legend(axT, [labels, {'threshold'}], 'Location', 'southeast', 'AutoUpdate', 'off', ...
    'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);

% Bottom: three stat bars.
bar_panel(subplot(2, 3, 4), labels, avail,    col_of, 'availability [%]',  is_oracle, false);
bar_panel(subplot(2, 3, 5), labels, deadtime, col_of, 'dead time [s]',     is_oracle, false);
bar_panel(subplot(2, 3, 6), labels, gap,      col_of, 'mean oracle gap [dB]', is_oracle, true);

exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function shade_jammer_on(ax, scn)
% Light-gray vertical bands over the intervals where the jammer transmits.
on = scn.jammer_on(:)';
d  = diff([0, on, 0]);
starts = find(d == 1);
stops  = find(d == -1) - 1;
yl = [-1e3, 1e3];
for i = 1:numel(starts)
    x0 = scn.t_s(starts(i)); x1 = scn.t_s(stops(i));
    patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end


function bar_panel(ax, labels, vals, col_of, ylab, is_oracle, drop_oracle)
% One labelled bar chart, bars colored to match the timeline legend.
hold(ax, 'on'); grid(ax, 'on');
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
idx = 1:numel(labels);
if drop_oracle, idx = find(~is_oracle); end
for j = 1:numel(idx)
    i = idx(j);
    bar(ax, j, vals(i), 0.6, 'FaceColor', col_of(labels{i}), 'EdgeColor', 'none');
    if ~isnan(vals(i))
        text(ax, j, vals(i), sprintf('%.2f', vals(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8);
    end
end
set(ax, 'XTick', 1:numel(idx), 'XTickLabel', labels(idx));
ylabel(ax, ylab);
end


function exportgraphics_compat(fig, out_path)
% exportgraphics exists from R2020a; fall back to print on older releases.
if exist('exportgraphics', 'file')
    exportgraphics(fig, out_path, 'Resolution', 150);
else
    print(fig, out_path, '-dpng', '-r150');
end
end
