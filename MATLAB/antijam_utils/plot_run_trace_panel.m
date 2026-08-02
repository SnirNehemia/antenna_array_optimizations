function plot_run_trace_panel(ax_tr, scenario, run_log, aj, dir_s_db, k)
% PLOT_RUN_TRACE_PANEL  Draw the live SINR/directivity trace panel up to step k.
%
%   PLOT_RUN_TRACE_PANEL(ax_tr, scenario, run_log, aj, dir_s_db, k)
%
%   Achieved SINR, oracle SINR, the sinr_min_db threshold, and (right axis)
%   the directivity toward theta_s, plotted from step 1 to k. Shared by
%   save_run_gif (called once per animation frame) and save_run_trace_png
%   (called once, with k = end) so both stay visually identical.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

cla(ax_tr);
yyaxis(ax_tr, 'left'); cla(ax_tr); hold(ax_tr, 'on'); grid(ax_tr, 'on');
plot(ax_tr, scenario.t_s(1:k), run_log.sinr_db(1:k), 'b-', 'LineWidth', 1.2);
plot(ax_tr, scenario.t_s(1:k), run_log.oracle_sinr_db(1:k), 'k:', 'LineWidth', 1.0);
yline(ax_tr, aj.sinr_min_db, 'r--');
plot(ax_tr, scenario.t_s(k), run_log.sinr_db(k), 'bo', 'MarkerFaceColor', 'b');
ylabel(ax_tr, 'SINR [dB]');
ylim(ax_tr, [0, max(run_log.oracle_sinr_db)]);
yyaxis(ax_tr, 'right');
plot(ax_tr, scenario.t_s(1:k), dir_s_db(1:k), '-', 'Color', [0.9 0.5 0.1], ...
    'LineWidth', 1.0);
ylabel(ax_tr, 'directivity @ \theta_s [dBi]');
ylim(ax_tr, [min(dir_s_db) - 1, max(dir_s_db) + 1]);
xlim(ax_tr, [0, scenario.t_s(end)]);
xlabel(ax_tr, 't [s]');
legend(ax_tr, {'SINR', 'oracle', 'threshold'}, 'Location', 'southeast', ...
    'AutoUpdate', 'off', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);
end
