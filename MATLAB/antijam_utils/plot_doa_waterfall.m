function plot_doa_waterfall(predict_log, scn, aj, out_path)
% PLOT_DOA_WATERFALL  Diagnostics for the P8 predictive Mode C algorithm.
%
%   PLOT_DOA_WATERFALL(predict_log, scn, aj, out_path)
%
%   Makes the MUSIC + on/off-prediction mechanism visible. Four panels:
%     (1) theta tracking : true jammer theta_j(t) vs the MUSIC estimate
%         theta_hat(t) (drawn only where the jammer is detected present).
%     (2) phi tracking   : same for azimuth phi.
%     (3) presence       : true jammer-on(t) vs MUSIC-detected present(t) vs the
%         predicted-on(t) schedule that drives the pre-null.
%     (4) periodogram    : |FFT| of the presence record vs candidate period —
%         the "spectrogram of repeating patterns" the customer asked for — with
%         the detected toggle period and the true period marked.
%
%   NOTE: geometry is native 2-D (theta, phi) since [P7], so the DoA is shown as
%   two tracking traces rather than a single-angle pseudospectrum waterfall.
%
%   Inputs:
%       predict_log : run_log from closed_loop_run('predict', ...). Uses fields
%                     doa_theta_deg, doa_phi_deg, present, predicted_on,
%                     period_est (1 x T each), and sinr_db.
%       scn         : scenario struct from sim_scenario (truth + timing).
%       aj          : antijam config (unused today; kept for signature symmetry).
%       out_path    : PNG destination path.
%
%   Outputs: none (figure file on disk).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

t   = scn.t_s;
dt  = t(2) - t(1);
on  = scn.jammer_on(:)';
det = predict_log.present(:)';

fig = figure('Visible', 'off', 'Position', [50 50 1000 720], 'Color', 'w');

% ── (1) theta tracking ─────────────────────────────────────────────
ax1 = subplot(2, 2, 1); hold(ax1, 'on'); grid(ax1, 'on');
set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
th_true = scn.theta_j_deg(:)'; th_true(~on) = NaN;      % show truth only while ON
th_est  = predict_log.doa_theta_deg(:)';                 % already NaN when absent
plot(ax1, t, th_true, 'k-',  'LineWidth', 1.6);
plot(ax1, t, th_est,  'r.',  'MarkerSize', 6);
ylabel(ax1, '\theta_j [deg]'); xlabel(ax1, 't [s]');
title(ax1, 'MUSIC \theta estimate vs truth', 'Color', 'k');
legend(ax1, {'true', 'MUSIC'}, 'Location', 'best', 'AutoUpdate', 'off', ...
    'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);

% ── (2) phi tracking ───────────────────────────────────────────────
ax2 = subplot(2, 2, 2); hold(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
ph_true = scn.phi_j_deg(:)'; ph_true(~on) = NaN;
ph_est  = predict_log.doa_phi_deg(:)';
plot(ax2, t, ph_true, 'k-', 'LineWidth', 1.6);
plot(ax2, t, ph_est,  'r.', 'MarkerSize', 6);
ylabel(ax2, '\phi_j [deg]'); xlabel(ax2, 't [s]');
title(ax2, 'MUSIC \phi estimate vs truth', 'Color', 'k');
legend(ax2, {'true', 'MUSIC'}, 'Location', 'best', 'AutoUpdate', 'off', ...
    'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);

% ── (3) presence / prediction timeline ─────────────────────────────
ax3 = subplot(2, 2, 3); hold(ax3, 'on'); grid(ax3, 'on');
set(ax3, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
stairs(ax3, t, on,  'k-', 'LineWidth', 1.6);
stairs(ax3, t, det + 1.2, 'b-', 'LineWidth', 1.2);                 % offset for clarity
stairs(ax3, t, double(predict_log.predicted_on(:)') + 2.4, 'r-', 'LineWidth', 1.2);
set(ax3, 'YTick', [0.5 1.7 2.9], 'YTickLabel', {'true on', 'MUSIC present', 'predicted on'});
ylim(ax3, [-0.2 3.6]); xlabel(ax3, 't [s]');
title(ax3, 'Jammer presence: truth vs detected vs predicted', 'Color', 'k');

% ── (4) presence periodogram ───────────────────────────────────────
ax4 = subplot(2, 2, 4); hold(ax4, 'on'); grid(ax4, 'on');
set(ax4, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
p = det - mean(det);
Np = numel(p);
if any(p ~= 0)
    P = abs(fft(p)).^2;
    half = floor(Np / 2);
    bins = 2:half + 1;                                  % skip DC
    period_s = (Np ./ (bins - 1)) * dt;                 % bin -> period [s]
    plot(ax4, period_s, P(bins), 'b-', 'LineWidth', 1.2);
end
det_period_s = predict_log.period_est(end) * dt;        % learned period
if ~isnan(det_period_s)
    xline(ax4, det_period_s, 'r--', 'LineWidth', 1.4);
end
true_period_s = true_toggle_period(on, dt);
if ~isnan(true_period_s)
    xline(ax4, true_period_s, 'k:', 'LineWidth', 1.2);
end
xlabel(ax4, 'candidate period [s]'); ylabel(ax4, 'presence power');
title(ax4, 'Presence periodogram (on/off detection)', 'Color', 'k');
if ~isnan(det_period_s)
    xlim(ax4, [0, max(3 * det_period_s, 2 * dt)]);
end
legend(ax4, {'periodogram', 'detected', 'true'}, 'Location', 'best', ...
    'AutoUpdate', 'off', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);

exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function period_s = true_toggle_period(on, dt)
% Median spacing between jammer turn-on transitions (NaN if fewer than two).
starts = find(diff([0, on]) == 1);
if numel(starts) < 2
    period_s = NaN;
else
    period_s = median(diff(starts)) * dt;
end
end


function exportgraphics_compat(fig, out_path)
if exist('exportgraphics', 'file')
    exportgraphics(fig, out_path, 'Resolution', 150);
else
    print(fig, out_path, '-dpng', '-r150');
end
end
