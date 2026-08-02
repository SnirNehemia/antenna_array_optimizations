function save_run_trace_png(run_log, scenario, stack1, stack2, theta_deg, phi_deg, ...
                            aj, title_label, out_path)
% SAVE_RUN_TRACE_PNG  Fast static PNG of one closed-loop run's live SINR trace.
%
%   SAVE_RUN_TRACE_PNG(run_log, scenario, stack1, stack2, theta_deg, phi_deg, ...
%                      aj, title_label, out_path)
%
%   Draws only the bottom trace panel of SAVE_RUN_GIF (achieved SINR, oracle
%   SINR, the sinr_min_db threshold, and directivity toward theta_s), fully
%   populated for the whole run, and saves it as one PNG. There is no
%   per-frame loop, no pattern heatmap, and no video encoding — video
%   rendering (SAVE_RUN_GIF) dominates the demo's runtime because it redraws
%   and re-encodes a pcolor heatmap for up to gif_cfg.max_frames frames, while
%   this trace panel is just line plots and costs one figure render. Use this
%   for a near-instant "what happened" glimpse (dead-time dips, threshold
%   crossings) without waiting on the full animated video.
%
%   Inputs: same as SAVE_RUN_GIF minus gif_cfg (no video config needed).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

dir_s_db = compute_directivity_trace(run_log, stack1, stack2, theta_deg, phi_deg);
T = numel(scenario.t_s);

fig = figure('Visible', 'off', 'Position', [50 50 760 280], 'Color', 'w');
ax_tr = axes(fig);
set(ax_tr, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
plot_run_trace_panel(ax_tr, scenario, run_log, aj, dir_s_db, T);
title(ax_tr, strrep(title_label, '_', '\_'), 'Interpreter', 'tex', 'Color', 'k');

exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function exportgraphics_compat(fig, out_path)
% exportgraphics exists from R2020a; fall back to print on older releases.
if exist('exportgraphics', 'file')
    exportgraphics(fig, out_path, 'Resolution', 150, 'BackgroundColor', 'w');
else
    print(fig, out_path, '-dpng', '-r150');
end
end
