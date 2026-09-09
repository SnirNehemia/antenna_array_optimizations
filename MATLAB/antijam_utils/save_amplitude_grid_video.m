function save_amplitude_grid_video(cases, stack1, stack2, theta_deg, phi_deg, ...
                                   aj, vid_cfg, title_label, out_path)
% SAVE_AMPLITUDE_GRID_VIDEO  One algorithm across four signal/jammer regimes.
%
%   SAVE_AMPLITUDE_GRID_VIDEO(cases, stack1, stack2, theta_deg, phi_deg, ...
%                             aj, vid_cfg, title_label, out_path)
%
%   Renders a 2x2 grid of (weak/strong desired signal) x (weak/strong jammer).
%   Each quadrant shows that regime's radiation pattern above its own output-SINR
%   trace, so the beam and what it bought sit together.
%
%   TWO SHARED SCALES, and they are the whole point of the figure:
%     * one colour scale across all four patterns, so a beam that has thrown away
%       gain looks different from one that has not;
%     * one SINR axis range across all four traces, so 25 dB in the strong-signal
%       quadrant is visibly 25 dB and not renormalised to look like the 5 dB in
%       the weak one. Per-panel autoscaling would make every regime look equally
%       healthy, which is precisely the impression this figure exists to prevent.
%
%   Inputs:
%       cases     : (1 x 4) struct array, read in this order ->
%                     (1,1) weak signal / weak jammer   (1,2) weak / strong
%                     (2,1) strong / weak               (2,2) strong / strong
%                   Fields per element:
%                     log        : struct from closed_loop_run for the algorithm
%                                  under test.
%                     oracle     : struct from closed_loop_run('oracle', ...).
%                     scn        : struct from sim_scenario. All four must share
%                                  the same time base.
%                     sigma_s_db : desired-signal power [dB re the noise floor].
%                     jn_ratio_db: jammer-to-noise ratio [dB].
%       stack1/2  : (N_el x N_theta x N_phi) far-field stacks (stack2 may be []).
%       theta_deg, phi_deg : full grids [deg].
%       aj        : antijam config (theta_s_deg, phi_s_deg, sinr_min_db).
%       vid_cfg   : struct. Required: max_frames, fps, dynamic_range_db.
%                   Optional: format 'mp4' (default) | 'gif'.
%       title_label : char naming the array and geometry.
%       out_path  : destination path.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [O2].

for key = {'max_frames', 'fps', 'dynamic_range_db'}
    if ~isfield(vid_cfg, key{1}) || isempty(vid_cfg.(key{1}))
        error('save_amplitude_grid_video:MissingKey', ...
            'Missing required video config key: ''%s''.', key{1});
    end
end
if numel(cases) ~= 4
    error('save_amplitude_grid_video:BadCases', ...
        'Expected 4 amplitude regimes; got %d.', numel(cases));
end
T = numel(cases(1).scn.t_s);
for i = 2:4
    if numel(cases(i).scn.t_s) ~= T
        error('save_amplitude_grid_video:TimeBase', ...
            ['All four regimes must share a time base so the frames line up; ' ...
             'regime %d has %d steps against %d.'], i, numel(cases(i).scn.t_s), T);
    end
end

fmt = 'mp4';
if isfield(vid_cfg, 'format') && ~isempty(vid_cfg.format), fmt = lower(vid_cfg.format); end

n_el    = size(stack1, 1);
n_theta = numel(theta_deg);
n_phi   = numel(phi_deg);
flat1   = reshape(stack1, n_el, []);
flat2   = [];
if ~isempty(stack2), flat2 = reshape(stack2, n_el, []); end

stride = max(1, ceil(T / vid_cfg.max_frames));
frames = 1:stride:T;

% ── Display-only decimation of the pattern grid ───────────────────
% Eight pcolor panels of a 181x360 grid is ~500k patches per frame, and
% getframe ran out of memory building the offscreen framebuffer for it. The
% panels are a few hundred pixels wide, so drawing every sample was never
% visible anyway: decimate to at most ~120x180 for DISPLAY. The underlying
% directivity is still computed on the full grid, so nothing quantitative
% changes -- only how many patches the renderer is asked to draw.
d_th = max(1, ceil(n_theta / 120));
d_ph = max(1, ceil(n_phi / 180));
i_th = 1:d_th:n_theta;
i_ph = 1:d_ph:n_phi;
theta_disp = theta_deg(i_th);
phi_disp   = phi_deg(i_ph);

% ── Shared scales ─────────────────────────────────────────────────
k_ref = frames(max(1, round(numel(frames) / 2)));
cmax = -Inf;
for i = 1:4
    d = frame_dbi(cases(i).log.W(:, k_ref), flat1, flat2, theta_deg, n_theta, n_phi);
    cmax = max(cmax, max(d(:)));
end
cmax = ceil(cmax) + 1;
dr   = vid_cfg.dynamic_range_db;

all_sinr = [];
for i = 1:4
    all_sinr = [all_sinr, cases(i).log.sinr_db(:)', cases(i).oracle.sinr_db(:)']; %#ok<AGROW>
end
y_lo = floor(min(all_sinr) / 5) * 5 - 2;
% Headroom for the legend strip: placed inside the axes it would otherwise sit
% on top of the traces, and the panels are too short to spare vertical space
% for an outside placement.
y_hi = ceil(max(all_sinr) / 5) * 5 + 2;
y_hi = y_hi + 0.22 * (y_hi - y_lo);

writer = [];
if strcmp(fmt, 'mp4')
    writer = VideoWriter(out_path, 'MPEG-4');
    writer.FrameRate = vid_cfg.fps;
    open(writer);
end

% Modest canvas on purpose: getframe allocates an offscreen framebuffer for
% the whole figure, and this machine runs several MATLAB sessions at once.
fig = figure('Visible', 'off', 'Position', [30 30 1040 860], 'Color', 'w');
% Rows 1 and 3 hold patterns, rows 2 and 4 the SINR traces beneath them, so each
% quadrant reads as one block.
ax_map = gobjects(1, 4);
ax_tr  = gobjects(1, 4);
slot_map = [1 2 5 6];
slot_tr  = [3 4 7 8];
for i = 1:4
    ax_map(i) = subplot(4, 2, slot_map(i));
    ax_tr(i)  = subplot(4, 2, slot_tr(i));
end
% Pull every axes down so the two header lines above have their own space, and
% squeeze each vertically so a pattern's phi label clears the panel beneath it.
for a = [ax_map, ax_tr]
    q = get(a, 'Position');
    set(a, 'Position', [q(1), q(2) * 0.93, q(3), q(4) * 0.82]);
end
set([ax_map, ax_tr], 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
    'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.15);

first = true;
fi = 0;
for k = frames
    fi = fi + 1;
    for i = 1:4
        c = cases(i);
        % ── pattern ───────────────────────────────────────────────
        cla(ax_map(i)); hold(ax_map(i), 'on');
        dbi = frame_dbi(c.log.W(:, k), flat1, flat2, theta_deg, n_theta, n_phi);
        pcolor(ax_map(i), phi_disp, theta_disp, dbi(i_th, i_ph));
        shading(ax_map(i), 'flat');
        colormap(ax_map(i), 'jet');
        caxis(ax_map(i), [cmax - dr, cmax]);
        set(ax_map(i), 'YDir', 'reverse');
        xlim(ax_map(i), [min(phi_deg), max(phi_deg)]);
        ylim(ax_map(i), [min(theta_deg), max(theta_deg)]);
        plot(ax_map(i), aj.phi_s_deg, aj.theta_s_deg, 'gp', ...
            'MarkerSize', 12, 'MarkerFaceColor', 'g');
        if c.scn.jammer_on(k) && mod(fi, 2) == 1
            scatter(ax_map(i), c.scn.phi_j_deg(k), c.scn.theta_j_deg(k), 170, 'o', ...
                'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r', ...
                'MarkerEdgeAlpha', 0.6, 'LineWidth', 2);
        end
        title(ax_map(i), sprintf('\\bf%s signal \\sigma_s = %.0f dB   \\bf%s jammer J/N = %.0f dB', ...
            strength_word(c.sigma_s_db, 10), c.sigma_s_db, ...
            strength_word(c.jn_ratio_db, 20), c.jn_ratio_db), ...
            'Interpreter', 'tex', 'Color', 'k');
        ylabel(ax_map(i), '\theta [deg]');
        % No xlabel on the pattern panels. It sat immediately above the SINR
        % panel's own title and the two overlapped. The axis meaning is stated
        % once in the header line instead, which cannot collide with anything.
        if mod(i, 2) == 0
            cb = colorbar(ax_map(i));
            cb.Label.String = 'directivity [dBi]';
            set(cb, 'Color', 'k');
            set(cb.Label, 'Color', 'k');
        end

        % ── SINR trace, shared range ──────────────────────────────
        cla(ax_tr(i)); hold(ax_tr(i), 'on');
        shade_on_windows(ax_tr(i), c.scn, y_lo, y_hi);
        plot(ax_tr(i), c.scn.t_s(1:k), c.oracle.sinr_db(1:k), '--', ...
            'Color', [0.45 0.45 0.45], 'LineWidth', 1.8);
        plot(ax_tr(i), c.scn.t_s(1:k), c.log.sinr_db(1:k), '-', ...
            'Color', [0.106 0.686 0.478], 'LineWidth', 2);
        yline(ax_tr(i), aj.sinr_min_db, ':', 'Color', [0.75 0.2 0.2], 'LineWidth', 1.4);
        xlim(ax_tr(i), [c.scn.t_s(1), c.scn.t_s(end)]);
        ylim(ax_tr(i), [y_lo, y_hi]);          % SHARED across all four
        grid(ax_tr(i), 'on');
        ylabel(ax_tr(i), 'SINR [dB]');
        if i >= 3, xlabel(ax_tr(i), 'time [s]'); end
        % The achieved/achievable pair as the SINR panel's own title. An earlier
        % attempt to place it inside the axes with text(ax,'Units','normalized',
        % 'Position',...) landed in the wrong panel entirely, so it lives here,
        % where placement is unambiguous -- the collision it used to have with
        % the pattern's phi label was removed by dropping that label instead.
        title(ax_tr(i), sprintf('\\rm\\fontsize{9}achieved %.1f dB   ·   achievable %.1f dB', ...
            c.log.sinr_db(k), c.oracle.sinr_db(k)), 'Interpreter', 'tex', 'Color', 'k');
        if i == 1
            lg = legend(ax_tr(i), {'achievable (oracle)', 'achieved', 'threshold'}, ...
                'Orientation', 'horizontal', 'Location', 'northwest', ...
                'Box', 'off', 'Interpreter', 'none');
            set(lg, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'none', 'FontSize', 8);
        end
    end

    if cases(1).scn.jammer_on(k)
        st = 'jammer ON';
    else
        st = 'jammer OFF';
    end
    % An annotation at an explicit normalized position, not sgtitle: sgtitle
    % places itself above the axes and was clipped against the figure edge.
    delete(findall(fig, 'Type', 'annotation'));
    annotation(fig, 'textbox', [0.02 0.955 0.96 0.042], ...
        'String', sprintf('%s  —  t = %.1f s  —  %s', ...
            strrep(title_label, '_', '\_'), cases(1).scn.t_s(k), st), ...
        'Interpreter', 'tex', 'Color', 'k', 'FontWeight', 'bold', ...
        'FontSize', 12, 'HorizontalAlignment', 'center', ...
        'LineStyle', 'none', 'VerticalAlignment', 'middle');
    annotation(fig, 'textbox', [0.02 0.925 0.96 0.030], ...
        'String', ['pattern axes: x = azimuth phi, y = elevation theta (deg)  ·  ' ...
                   'all powers re the per-element noise floor (0 dB)  ·  ' ...
                   'SINR axes share one range  ·  shaded = jammer ON'], ...
        'Interpreter', 'none', 'Color', 'k', 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'LineStyle', 'none', ...
        'VerticalAlignment', 'middle');

    set(findall(fig, 'Type', 'text'), 'Color', 'k');
    im = frame2im(getframe(fig));
    if strcmp(fmt, 'mp4')
        im = im(1:end - mod(size(im, 1), 2), 1:end - mod(size(im, 2), 2), :);
        writeVideo(writer, im);
    else
        [A, map] = rgb2ind(im, 256);
        if first
            imwrite(A, map, out_path, 'gif', 'LoopCount', Inf, 'DelayTime', 1 / vid_cfg.fps);
            first = false;
        else
            imwrite(A, map, out_path, 'gif', 'WriteMode', 'append', 'DelayTime', 1 / vid_cfg.fps);
        end
    end
end
if ~isempty(writer), close(writer); end
close(fig);
end


% NOTE ON `caxis`: kept deliberately. `clim` is the modern spelling but arrived in
% R2022a, and this project targets R2020a.
% ────────────────────────── HELPERS ───────────────────────────────

function w = strength_word(value_db, pivot_db)
if value_db >= pivot_db
    w = 'STRONG';
else
    w = 'weak';
end
end


function shade_on_windows(ax, scenario, y_lo, y_hi)
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


function dbi = frame_dbi(w, flat1, flat2, theta_deg, n_theta, n_phi)
p = abs(w' * flat1).^2;
if ~isempty(flat2)
    p = p + abs(w' * flat2).^2;
end
p = reshape(p, n_theta, n_phi);
sw = sind(theta_deg(:)) * ones(1, n_phi);
p_avg = sum(p(:) .* sw(:)) / sum(sw(:));
dbi = 10 * log10(max(p, realmin) / max(p_avg, realmin));
end
