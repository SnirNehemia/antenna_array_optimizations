function save_comparison_video(run_logs, labels, scenario, stack1, stack2, ...
                               theta_deg, phi_deg, aj, vid_cfg, title_label, out_path)
% SAVE_COMPARISON_VIDEO  Several algorithms on one scenario, side by side.
%
%   SAVE_COMPARISON_VIDEO(run_logs, labels, scenario, stack1, stack2, ...
%                         theta_deg, phi_deg, aj, vid_cfg, title_label, out_path)
%
%   `save_run_gif` animates ONE run, which cannot answer "why is this algorithm
%   better than that one" -- the viewer has to remember the other file. This
%   renders every algorithm's radiation pattern on a shared colour scale in one
%   row, above a single SINR panel carrying all of their traces plus the oracle
%   and the operating threshold. Differences read as differences in the picture
%   rather than between two files.
%
%   The shared colour scale is the point and not a convenience: per-panel
%   autoscaling would make a beamformer that has thrown away 10 dB of gain look
%   identical to one that has not.
%
%   Inputs:
%       run_logs    : (1 x M) cell of structs from closed_loop_run, each with
%                     oracle_sinr_db added. All must share the scenario. Do NOT
%                     include the oracle run itself -- it is drawn from
%                     run_logs{1}.oracle_sinr_db as the dashed reference, and
%                     passing it again puts it on the plot and in the legend
%                     twice.
%       labels      : (1 x M) cellstr, the algorithm name for each panel.
%       scenario    : struct from sim_scenario.
%       stack1/2    : (N_el x N_theta x N_phi) far-field stacks (stack2 may be []).
%       theta_deg, phi_deg : full grids [deg].
%       aj          : antijam config (theta_s_deg, phi_s_deg, sinr_min_db).
%       vid_cfg     : struct. Required: max_frames, fps, dynamic_range_db.
%                     Optional: format 'mp4' (default) | 'gif'.
%       title_label : char shown in the figure title (e.g. 'spacing0.6 / T=10 s').
%       out_path    : destination path; the extension should match the format.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [O].

for key = {'max_frames', 'fps', 'dynamic_range_db'}
    if ~isfield(vid_cfg, key{1}) || isempty(vid_cfg.(key{1}))
        error('save_comparison_video:MissingKey', ...
            'Missing required video config key: ''%s''.', key{1});
    end
end
if numel(run_logs) ~= numel(labels)
    error('save_comparison_video:LabelMismatch', ...
        'Got %d run logs but %d labels.', numel(run_logs), numel(labels));
end
fmt = 'mp4';
if isfield(vid_cfg, 'format') && ~isempty(vid_cfg.format), fmt = lower(vid_cfg.format); end

M       = numel(run_logs);
T       = numel(scenario.t_s);
stride  = max(1, ceil(T / vid_cfg.max_frames));
frames  = 1:stride:T;
n_el    = size(stack1, 1);
n_theta = numel(theta_deg);
n_phi   = numel(phi_deg);
flat1   = reshape(stack1, n_el, []);
flat2   = [];
if ~isempty(stack2), flat2 = reshape(stack2, n_el, []); end

% ── One colour scale for every panel and every frame ──────────────
% Taken as the max over ALL algorithms at a mid-run frame, so a panel that has
% lost gain shows it instead of being renormalised back to full scale.
k_ref = frames(max(1, round(numel(frames) / 2)));
cmax = -Inf;
for m = 1:M
    d = frame_dbi(run_logs{m}.W(:, k_ref), flat1, flat2, theta_deg, phi_deg, ...
        n_theta, n_phi);
    cmax = max(cmax, max(d(:)));
end
cmax = ceil(cmax) + 1;
dr   = vid_cfg.dynamic_range_db;

% SINR traces, plus the oracle from the first log (identical across logs).
sinr = zeros(M, T);
for m = 1:M, sinr(m, :) = run_logs{m}.sinr_db(:)'; end
o_sinr = run_logs{1}.oracle_sinr_db(:)';
y_lo = min([sinr(:); o_sinr(:)]) - 2;
y_hi = max([sinr(:); o_sinr(:)]) + 2;

% Categorical colours, fixed per algorithm (never by rank).
series = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; ...
          0.929 0.631 0.000];

writer = [];
if strcmp(fmt, 'mp4')
    writer = VideoWriter(out_path, 'MPEG-4');
    writer.FrameRate = vid_cfg.fps;
    open(writer);
end

W_fig = min(420 * M + 120, 1700);
fig = figure('Visible', 'off', 'Position', [40 40 W_fig 640], 'Color', 'w');
ax_map = gobjects(1, M);
for m = 1:M
    ax_map(m) = subplot(2, M, m);
end
ax_tr = subplot(2, 1, 2);
set([ax_map, ax_tr], 'Color', 'w', 'XColor', 'k', 'YColor', 'k');

first = true;
fi = 0;
for k = frames
    fi = fi + 1;
    th_j = scenario.theta_j_deg(k);
    ph_j = scenario.phi_j_deg(k);
    jam_on = scenario.jammer_on(k);

    for m = 1:M
        cla(ax_map(m)); hold(ax_map(m), 'on');
        dbi = frame_dbi(run_logs{m}.W(:, k), flat1, flat2, theta_deg, phi_deg, ...
            n_theta, n_phi);
        pcolor(ax_map(m), phi_deg, theta_deg, dbi);
        shading(ax_map(m), 'flat');
        colormap(ax_map(m), 'jet');
        caxis(ax_map(m), [cmax - dr, cmax]);
        set(ax_map(m), 'YDir', 'reverse');
        xlim(ax_map(m), [min(phi_deg), max(phi_deg)]);
        ylim(ax_map(m), [min(theta_deg), max(theta_deg)]);
        plot(ax_map(m), aj.phi_s_deg, aj.theta_s_deg, 'gp', ...
            'MarkerSize', 13, 'MarkerFaceColor', 'g');
        if jam_on && mod(fi, 2) == 1
            scatter(ax_map(m), ph_j, th_j, 190, 'o', 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', 'r', 'MarkerEdgeAlpha', 0.6, 'LineWidth', 2);
        end
        % The achieved SINR in each panel's own title: the picture shows the
        % pattern, the number says what it bought.
        title(ax_map(m), sprintf('%s  (%.1f dB)', ...
            strrep(labels{m}, '_', '\_'), sinr(m, k)), ...
            'Interpreter', 'tex', 'Color', 'k', 'FontWeight', 'bold');
        if m == 1
            ylabel(ax_map(m), 'Elevation \theta [deg]');
        end
        xlabel(ax_map(m), 'Azimuth \phi [deg]');
        if m == M
            cb = colorbar(ax_map(m));
            cb.Label.String = 'directivity [dBi]';
        end
    end

    % ── Shared SINR panel ─────────────────────────────────────────
    cla(ax_tr); hold(ax_tr, 'on');
    % Shade the ON windows so the toggling is legible without reading the title.
    shade_on_windows(ax_tr, scenario, y_lo, y_hi);
    plot(ax_tr, scenario.t_s(1:k), o_sinr(1:k), '--', 'Color', [0.45 0.45 0.45], ...
        'LineWidth', 2);
    for m = 1:M
        plot(ax_tr, scenario.t_s(1:k), sinr(m, 1:k), '-', ...
            'Color', series(mod(m - 1, size(series, 1)) + 1, :), 'LineWidth', 2);
    end
    yline(ax_tr, aj.sinr_min_db, ':', 'Color', [0.75 0.2 0.2], 'LineWidth', 1.5);
    xlim(ax_tr, [scenario.t_s(1), scenario.t_s(end)]);
    ylim(ax_tr, [y_lo, y_hi]);
    xlabel(ax_tr, 'time [s]'); ylabel(ax_tr, 'output SINR [dB]');
    legend(ax_tr, [{'oracle'}, labels(:)', {'threshold'}], ...
        'Location', 'southeast', 'Box', 'off', 'Interpreter', 'none');
    grid(ax_tr, 'on');
    if jam_on
        st = sprintf('jammer ON (J/N %.0f dB)', scenario.jn_ratio_db(k));
    else
        st = 'jammer OFF';
    end
    title(ax_tr, sprintf('%s — t = %.1f s — %s', ...
        strrep(title_label, '_', '\_'), scenario.t_s(k), st), ...
        'Interpreter', 'tex', 'Color', 'k');

    im = frame2im(getframe(fig));
    if strcmp(fmt, 'mp4')
        im = im(1:end - mod(size(im, 1), 2), 1:end - mod(size(im, 2), 2), :);
        writeVideo(writer, im);
    else
        [A, map] = rgb2ind(im, 256);
        if first
            imwrite(A, map, out_path, 'gif', 'LoopCount', Inf, ...
                'DelayTime', 1 / vid_cfg.fps);
            first = false;
        else
            imwrite(A, map, out_path, 'gif', 'WriteMode', 'append', ...
                'DelayTime', 1 / vid_cfg.fps);
        end
    end
end
if ~isempty(writer), close(writer); end
close(fig);
end


% NOTE ON `caxis`: kept deliberately. `clim` is the modern spelling but was
% introduced in R2022a, and this project targets R2020a (see docs/notes.md), so
% the Code Analyzer's suggestion would break the supported baseline.
% ────────────────────────── HELPERS ───────────────────────────────

function shade_on_windows(ax, scenario, y_lo, y_hi)
% Light vertical bands over the intervals where the jammer transmits.
on = scenario.jammer_on(:)';
d  = diff([false, on, false]);
starts = find(d == 1);
stops  = find(d == -1) - 1;
t = scenario.t_s;
for i = 1:numel(starts)
    x0 = t(starts(i));
    x1 = t(min(stops(i), numel(t)));
    patch(ax, [x0 x1 x1 x0], [y_lo y_lo y_hi y_hi], [0.93 0.93 0.93], ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end


function dbi = frame_dbi(w, flat1, flat2, theta_deg, ~, n_theta, n_phi)
% Directivity dBi grid for one weight vector (power-summed over components),
% normalised by the solid-angle-weighted mean, per project convention.
p = abs(w' * flat1).^2;
if ~isempty(flat2)
    p = p + abs(w' * flat2).^2;
end
p = reshape(p, n_theta, n_phi);
sw = sind(theta_deg(:)) * ones(1, n_phi);
p_avg = sum(p(:) .* sw(:)) / sum(sw(:));
dbi = 10 * log10(max(p, realmin) / max(p_avg, realmin));
end
