function plot_freq_waterfall(freq_log, scn, notch_cfg, out_path)
% PLOT_FREQ_WATERFALL  Diagnostics for the P10 carrier tracker + RF notch.
%
%   PLOT_FREQ_WATERFALL(freq_log, scn, notch_cfg, out_path)
%
%   The spectral counterpart of plot_doa_waterfall. Three panels:
%     (1) SPECTROGRAM — the incoherently averaged periodogram vs time, with the
%         true carrier f_j(t) overlaid, the tracker's estimate on top of it, and
%         the notch's 3-dB band shaded. This is the figure that answers the
%         customer's actual question: which frequency is the jammer on.
%     (2) ESTIMATION ERROR — f_hat - f_j in kHz, against the notch half-width so
%         the error can be read directly as "well inside the notch" or not.
%         Jammer-OFF stretches are shaded: there the tracker is deliberately
%         HOLDING its last estimate, so the error there is not scored.
%     (3) SINR — delivered (notch engaged) vs the spatial-null-only baseline
%         recorded alongside it. NOTE these two traces coincide whenever the
%         spatial null is converged: the null has already driven the jammer far
%         below the noise floor and a notch can only remove what the jammer
%         still contributes. The gap opens exactly where the null is weak or has
%         not yet formed — turn-on transients, carrier/angle hops, and jammers
%         too close to the main beam to null. That is the honest reading of this
%         panel, and it is the point of the phase.
%
%   Inputs:
%       freq_log  : run_log from closed_loop_run with the [P10] waveform layer
%                   on. Uses f_hat_hz, freq_present, sinr_db, sinr_no_notch_db,
%                   notch_f_hz, and (for panel 1) pspec / pspec_k / pspec_f_hz.
%       scn       : scenario struct from sim_scenario (truth: f_j_hz, jammer_on).
%       notch_cfg : the notch config section (bw_hz used for the shaded band);
%                   may be [] to omit the band.
%       out_path  : PNG destination path.
%
%   Outputs: none (figure file on disk).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P10].

t  = scn.t_s;
on = scn.jammer_on(:).';
if ~isfield(freq_log, 'f_hat_hz')
    error('plot_freq_waterfall:NoFreqLog', ...
        'freq_log has no f_hat_hz — the run was made without a [P10] carrier tracker.');
end
half_bw_hz = NaN;
if ~isempty(notch_cfg) && isfield(notch_cfg, 'bw_hz')
    half_bw_hz = notch_cfg.bw_hz / 2;
end

fig = figure('Visible', 'off', 'Position', [50 50 1000 860], 'Color', 'w');

% ── (1) spectrogram + true carrier + estimate + notch band ─────────
ax1 = subplot(3, 1, 1); hold(ax1, 'on');
set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
if isfield(freq_log, 'pspec') && ~isempty(freq_log.pspec)
    tk = t(freq_log.pspec_k);
    fk = freq_log.pspec_f_hz / 1e6;
    S  = 10 * log10(max(freq_log.pspec, realmin));
    imagesc(ax1, tk, fk, S);
    set(ax1, 'YDir', 'normal');
    colormap(ax1, 'parula');
    cb = colorbar(ax1); cb.Color = 'k';
    ylabel(cb, 'periodogram [dB]');
    % Stretch the colour scale over the top of the range so the line stands out
    % against the wideband desired signal rather than being washed into it.
    % A fixed 40 dB dynamic range below the peak (the same convention the demo
    % videos use), so the carrier line stands out instead of being washed into
    % the wideband desired signal's floor.
    % caxis, NOT clim: clim is R2022a+ and the milestone targets R2020a
    % (docs/MATLAB_R2020a_changes.md). The Code Analyzer's "use clim" hint is
    % correct for modern MATLAB and wrong for this project — do not apply it.
    caxis(ax1, [max(S(:)) - 40, max(S(:))]); %#ok<CAXIS>
end
% Estimate first, truth on top — truth is the reference and must stay readable
% where the two coincide (which is almost everywhere, and is the point).
plot(ax1, t, freq_log.f_hat_hz / 1e6, 'r.', 'MarkerSize', 4);
f_true = scn.f_j_hz(:).' / 1e6; f_true(~on) = NaN;    % truth only while ON
plot(ax1, t, f_true, 'w--', 'LineWidth', 1.8);
if ~isnan(half_bw_hz)
    plot(ax1, t, (freq_log.notch_f_hz - half_bw_hz) / 1e6, 'r-', 'LineWidth', 0.5);
    plot(ax1, t, (freq_log.notch_f_hz + half_bw_hz) / 1e6, 'r-', 'LineWidth', 0.5);
end
xlim(ax1, [t(1), t(end)]);
ylim(ax1, [-1, 1] * max(abs(scn.f_j_hz)) / 1e6 * 1.6);
xlabel(ax1, 't [s]'); ylabel(ax1, 'offset from f_c [MHz]');
title(ax1, 'Jammer carrier: periodogram, truth, and tracked estimate', 'Color', 'k');
legend(ax1, {'estimate', 'true f_j'}, 'Location', 'northeast', 'AutoUpdate', 'off', ...
    'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);

% ── (2) estimation error vs the notch half-width ───────────────────
ax2 = subplot(3, 1, 2); hold(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
err_khz = (freq_log.f_hat_hz - scn.f_j_hz) / 1e3;
err_khz(~on) = NaN;                                   % not scored while silent
plot(ax2, t, err_khz, 'b-', 'LineWidth', 1.0);
if ~isnan(half_bw_hz)
    % The notch half-width is the scale that matters: error inside these lines
    % costs essentially none of the notch's rejection.
    yline_compat(ax2,  half_bw_hz / 1e3, 'r--');
    yline_compat(ax2, -half_bw_hz / 1e3, 'r--');
    ylim(ax2, [-1.4, 1.4] * half_bw_hz / 1e3);
end
xlim(ax2, [t(1), t(end)]);
shade_off(ax2, t, on);
xlabel(ax2, 't [s]'); ylabel(ax2, 'f_{hat} - f_j [kHz]');
title(ax2, 'Carrier-estimation error vs notch half-width (grey = jammer off, estimate held)', ...
    'Color', 'k');

% ── (3) SINR with the notch vs spatial null only ───────────────────
ax3 = subplot(3, 1, 3); hold(ax3, 'on'); grid(ax3, 'on');
set(ax3, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
if isfield(freq_log, 'sinr_no_notch_db')
    plot(ax3, t, freq_log.sinr_no_notch_db, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4);
end
plot(ax3, t, freq_log.sinr_db, 'b-', 'LineWidth', 1.2);
xlim(ax3, [t(1), t(end)]);
shade_off(ax3, t, on);
xlabel(ax3, 't [s]'); ylabel(ax3, 'SINR [dB]');
title(ax3, 'SINR: spatial null only vs null + RF notch', 'Color', 'k');
legend(ax3, {'spatial null only', 'null + notch'}, 'Location', 'best', ...
    'AutoUpdate', 'off', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);

exportgraphics_compat(fig, out_path);
close(fig);
end


% ────────────────────────── HELPERS ───────────────────────────────

function shade_off(ax, t, on)
% Light grey patches over the jammer-OFF stretches. MUST be called AFTER the
% data is plotted: the patches are sized to the axis limits the data produced,
% then pushed behind it, and the limits are restored (a patch spanning a huge
% y-range would otherwise blow the auto-scaling and flatten every trace).
yl = get(ax, 'YLim');
starts = find(diff([true, ~on]) == 1);
stops  = find(diff([~on, true]) == -1);
for i = 1:min(numel(starts), numel(stops))
    p = patch(ax, t([starts(i) stops(i) stops(i) starts(i)]), yl([1 1 2 2]), ...
        [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(p, 'bottom');
end
set(ax, 'YLim', yl);
end


function yline_compat(ax, y, spec)
% yline() is R2018b+, but keep the drawing style consistent with the repo's
% existing xline usage by falling back to a plain plot on the current x-limits.
xl = get(ax, 'XLim');
plot(ax, xl, [y y], spec, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end


function exportgraphics_compat(fig, out_path)
if exist('exportgraphics', 'file')
    exportgraphics(fig, out_path, 'Resolution', 150);
else
    print(fig, out_path, '-dpng', '-r150');
end
end
