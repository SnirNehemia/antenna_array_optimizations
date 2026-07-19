function save_run_gif(run_log, scenario, stack1, stack2, theta_deg, phi_deg, ...
                      aj, gif_cfg, title_label, out_path)
% SAVE_RUN_GIF  Animated GIF of one closed-loop run: 2-D pattern + live KPIs.
%
%   SAVE_RUN_GIF(run_log, scenario, stack1, stack2, theta_deg, phi_deg, ...
%                aj, gif_cfg, title_label, out_path)
%
%   Each frame shows:
%     - the full theta x phi radiation pattern (directivity dBi heatmap, same
%       pcolor/jet/reversed-theta style as the manual weight-tuner GUI) for
%       the weights applied at that instant, with a fixed color scale;
%     - the true jammer position as a dot on the map — filled magenta while
%       transmitting (marker size grows with its linear amplitude 10^(JNR/20)),
%       hollow gray while silent; the steered direction theta_s as a green
%       pentagram;
%     - a lower panel with the live traces up to the current instant:
%       achieved SINR, oracle SINR, the sinr_min_db threshold, and (right
%       axis) the noise-normalized gain toward theta_s.
%
%   Inputs:
%       run_log     : struct from closed_loop_run (+ oracle_sinr_db).
%       scenario    : struct from sim_scenario.
%       stack1/2    : full (N_el, N_theta, N_phi) stacks (stack2 may be []).
%       theta_deg, phi_deg : full grids [deg].
%       aj          : antijam config (theta_s_deg, sinr_min_db, cut_type,
%                     cut_theta_deg / cut_phi_deg).
%       gif_cfg     : struct. Required: max_frames, fps, dynamic_range_db.
%       title_label : char shown in the title (e.g. 'J1 / lcmv').
%       out_path    : destination .gif path.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

for key = {'max_frames', 'fps', 'dynamic_range_db'}
    if ~isfield(gif_cfg, key{1}) || isempty(gif_cfg.(key{1}))
        error('save_run_gif:MissingKey', 'Missing required gif config key: ''%s''.', key{1});
    end
end

T       = numel(scenario.t_s);
stride  = max(1, ceil(T / gif_cfg.max_frames));
frames  = 1:stride:T;
n_el    = size(stack1, 1);
n_theta = numel(theta_deg);
n_phi   = numel(phi_deg);
flat1   = reshape(stack1, n_el, []);
flat2   = [];
if ~isempty(stack2), flat2 = reshape(stack2, n_el, []); end

% Precompute the live traces (full resolution — cheap).
e_s = run_log.cut.e_s;
n_c = size(e_s, 2);
gain_db = zeros(1, T);
for k = 1:T
    w = run_log.W(:, k);
    gain_db(k) = 10 * log10((sum(abs(w' * e_s).^2) / n_c) / real(w' * w));
end

% Jammer / steer positions on the full grid.
[th_s, ph_s] = cut_to_theta_phi(aj.theta_s_deg, aj);

% Fixed color scale from the first frame's pattern.
dbi = frame_dbi(run_log.W(:, frames(1)), flat1, flat2, theta_deg, phi_deg, n_theta, n_phi);
cmax = ceil(max(dbi(:))) + 1;
dr   = gif_cfg.dynamic_range_db;

fig = figure('Visible', 'off', 'Position', [50 50 760 620], 'Color', 'w');
ax_map = subplot(3, 1, [1 2]);
ax_tr  = subplot(3, 1, 3);
% Force light styling regardless of the session's figure theme.
set([ax_map, ax_tr], 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
first = true;
for k = frames
    w = run_log.W(:, k);
    dbi = frame_dbi(w, flat1, flat2, theta_deg, phi_deg, n_theta, n_phi);

    % ── Pattern map ───────────────────────────────────────────────
    cla(ax_map); hold(ax_map, 'on');
    pcolor(ax_map, phi_deg, theta_deg, dbi); shading(ax_map, 'flat');
    colormap(ax_map, 'jet'); caxis(ax_map, [cmax - dr, cmax]);
    cb = colorbar(ax_map); cb.Label.String = 'directivity [dBi]';
    set(ax_map, 'YDir', 'reverse');
    xlim(ax_map, [min(phi_deg), max(phi_deg)]);
    ylim(ax_map, [min(theta_deg), max(theta_deg)]);
    xlabel(ax_map, 'Azimuth \phi [deg]'); ylabel(ax_map, 'Elevation \theta [deg]');
    plot(ax_map, ph_s, th_s, 'gp', 'MarkerSize', 14, 'MarkerFaceColor', 'g');
    [th_j, ph_j] = cut_to_theta_phi(scenario.theta_j_deg(k), aj);
    if scenario.jammer_on(k)
        msize = 6 + 0.7 * 10^(scenario.jn_ratio_db(k) / 20);
        plot(ax_map, ph_j, th_j, 'o', 'MarkerSize', msize, ...
            'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        state = sprintf('ON, J/N %.0f dB', scenario.jn_ratio_db(k));
    else
        plot(ax_map, ph_j, th_j, 'o', 'MarkerSize', 8, ...
            'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.5);
        state = 'OFF';
    end
    title(ax_map, sprintf('%s — t = %.1f s | jammer @ %.0f\\circ (%s)', ...
        strrep(title_label, '_', '\_'), scenario.t_s(k), ...
        scenario.theta_j_deg(k), state), 'Interpreter', 'tex');

    % ── Live traces ───────────────────────────────────────────────
    cla(ax_tr);
    yyaxis(ax_tr, 'left'); cla(ax_tr); hold(ax_tr, 'on'); grid(ax_tr, 'on');
    plot(ax_tr, scenario.t_s(1:k), run_log.sinr_db(1:k), 'b-', 'LineWidth', 1.2);
    plot(ax_tr, scenario.t_s(1:k), run_log.oracle_sinr_db(1:k), 'k:', 'LineWidth', 1.0);
    yline(ax_tr, aj.sinr_min_db, 'r--');
    plot(ax_tr, scenario.t_s(k), run_log.sinr_db(k), 'bo', 'MarkerFaceColor', 'b');
    ylabel(ax_tr, 'SINR [dB]');
    ylim(ax_tr, [min(-5, min(run_log.sinr_db) - 2), max(run_log.oracle_sinr_db) + 3]);
    yyaxis(ax_tr, 'right');
    plot(ax_tr, scenario.t_s(1:k), gain_db(1:k), '-', 'Color', [0.9 0.5 0.1], ...
        'LineWidth', 1.0);
    ylabel(ax_tr, 'gain @ \theta_s [dB]');
    ylim(ax_tr, [min(gain_db) - 1, max(gain_db) + 1]);
    xlim(ax_tr, [0, scenario.t_s(end)]);
    xlabel(ax_tr, 't [s]');
    legend(ax_tr, {'SINR', 'oracle', 'threshold'}, 'Location', 'southeast', ...
        'AutoUpdate', 'off');

    % ── Append frame ──────────────────────────────────────────────
    frame = getframe(fig);
    [A, map] = rgb2ind(frame2im(frame), 256);
    if first
        imwrite(A, map, out_path, 'gif', 'LoopCount', Inf, ...
            'DelayTime', 1 / gif_cfg.fps);
        first = false;
    else
        imwrite(A, map, out_path, 'gif', 'WriteMode', 'append', ...
            'DelayTime', 1 / gif_cfg.fps);
    end
end
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function dbi = frame_dbi(w, flat1, flat2, theta_deg, phi_deg, n_theta, n_phi)
% Directivity dBi grid for one weight vector (power-summed over components).
p = abs(w' * flat1).^2;
if ~isempty(flat2)
    p = p + abs(w' * flat2).^2;
end
af  = reshape(sqrt(p), n_theta, n_phi);   % magnitude AF; dbi squares it
dbi = compute_directivity_dbi_grid(af, theta_deg(:), phi_deg(:), []);
end


function [th, ph] = cut_to_theta_phi(angle_deg, aj)
% Map a cut-axis angle to (theta, phi) on the full grid.
if strcmp(aj.cut_type, 'phi_cut')
    th = aj.cut_theta_deg;
    ph = mod(angle_deg, 360);
else
    th = abs(angle_deg);
    ph = mod(aj.cut_phi_deg + 180 * (angle_deg < 0), 360);
end
end
