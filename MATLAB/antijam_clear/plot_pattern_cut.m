function figure_path = plot_pattern_cut(results, array, scenario, step_index, results_dir)
% ══════════════════════════════════════════════════════════════════
% PLOT_PATTERN_CUT
% The picture that makes the null visible.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   figure_path = PLOT_PATTERN_CUT(results, array, scenario, step_index, ...
%                                  results_dir)
%
%   One elevation cut of the radiation pattern, at the azimuth the signal and
%   jammer share, with one line per approach and three markings: the signal
%   direction, the jammer direction, and the guard sector.
%
%   The table says "12.1 dB". This says "a beam pointed at 30 degrees with a
%   notch at 55", which is a thing anyone in the room understands without being
%   told what a covariance is. It is worth its file for that alone.
%
%   THE ONE PLACE THE CONJUGATE APPEARS. This folder's convention is
%   y = w' * x, so the array gain in direction e is w' * e. matlab_utils/
%   compute_array_factor forms sum_n w_n * E_n, which is w.' * e. The two
%   differ by conj(w), and this is the single function that reconciles them --
%   see steering_vector for the full note. Everywhere else, w' * e holds.
%
%   Inputs:
%       results     : struct from run_closed_loop; uses names and final_weights.
%       array       : struct from make_array.
%       scenario    : struct from make_scenario. TRUTH -- for marking the true
%                     jammer direction on the figure. Reporting only.
%       step_index  : the step whose weights and jammer position are drawn.
%       results_dir : folder for the PNG. Created if absent.
%
%   Outputs:
%       figure_path : path of the saved PNG.

% ────────────────────────── CONSTANTS ─────────────────────────────

% Floor for the dB axis. Nulls go arbitrarily deep and an unbounded axis makes
% every other feature of the pattern invisible.
PATTERN_FLOOR_DBI = -40.0;

% Explicit colours everywhere. The MATLAB session theme overrides figure colour
% on print regardless of the property, so axes, text and legend are all set by
% hand rather than left to default.
LINE_COLOURS = [0.45 0.45 0.45;    % quiescent  - grey
                0.85 0.33 0.10;    % approach 1 - orange
                0.00 0.45 0.74];   % approach 2 - blue
FOREGROUND   = [0 0 0];
GUARD_SHADE  = [0.93 0.93 0.85];

theta_deg = array.theta_deg(:);
phi_index = nearest_index(array.phi_deg, array.profile.signal_phi_deg);

n_approaches = numel(results.names);

% ────────────────────────── PATTERN CUTS ──────────────────────────

pattern_dbi = nan(numel(theta_deg), n_approaches);
for i = 1:n_approaches
    weights = results.final_weights{i};

    % conj() converts this folder's w (for y = w'*x) into the convention
    % compute_array_factor expects (sum of w_n * E_n). See the header note.
    array_factor         = compute_array_factor(conj(weights), array.element_patterns);
    directivity_dbi_grid = compute_directivity_dbi_grid(array_factor, theta_deg, ...
                                                        array.phi_deg, []);
    pattern_dbi(:, i) = directivity_dbi_grid(:, phi_index);
end

pattern_dbi = max(pattern_dbi, PATTERN_FLOOR_DBI);

% ────────────────────────── DRAW ──────────────────────────────────

figure_handle = figure('Color', 'w', 'Position', [100 100 900 520]);
axes_handle   = axes('Parent', figure_handle, 'Color', 'w', ...
                     'XColor', FOREGROUND, 'YColor', FOREGROUND);
hold(axes_handle, 'on');

signal_theta_deg = scenario.signal_theta_deg;
jammer_theta_deg = scenario.jammer_theta_deg(step_index);
guard_deg        = array.profile.guard_deg;

% Guard sector: the region where a jammer is inside the main lobe and therefore
% out of scope. Drawn first so the pattern lines sit on top of it.
guard_left  = max(signal_theta_deg - guard_deg, theta_deg(1));
guard_right = min(signal_theta_deg + guard_deg, theta_deg(end));
patch('Parent', axes_handle, ...
      'XData', [guard_left guard_right guard_right guard_left], ...
      'YData', [PATTERN_FLOOR_DBI PATTERN_FLOOR_DBI 40 40], ...
      'FaceColor', GUARD_SHADE, 'EdgeColor', 'none');

plot_handles = gobjects(1, n_approaches);
for i = 1:n_approaches
    plot_handles(i) = plot(axes_handle, theta_deg, pattern_dbi(:, i), ...
                           'Color', LINE_COLOURS(i, :), 'LineWidth', 1.6);
end

y_limits = [PATTERN_FLOOR_DBI, max(pattern_dbi(:)) + 3];
plot(axes_handle, [signal_theta_deg signal_theta_deg], y_limits, ...
     'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.2);
plot(axes_handle, [jammer_theta_deg jammer_theta_deg], y_limits, ...
     'Color', [0.8 0 0], 'LineStyle', '--', 'LineWidth', 1.2);

% ────────────────────────── LABEL ─────────────────────────────────

legend_labels = [strrep(results.names, '_', ' '), {'signal', 'jammer'}];
legend_handle = legend(axes_handle, [plot_handles, ...
    plot(axes_handle, NaN, NaN, 'Color', [0 0.5 0], 'LineStyle', '--'), ...
    plot(axes_handle, NaN, NaN, 'Color', [0.8 0 0], 'LineStyle', '--')], ...
    legend_labels, 'Location', 'southwest');
set(legend_handle, 'TextColor', FOREGROUND, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

xlabel(axes_handle, 'elevation \theta [deg]', 'Color', FOREGROUND);
ylabel(axes_handle, 'directivity [dBi]', 'Color', FOREGROUND);
% Interpreter 'none': array names are plain text, not TeX.
title(axes_handle, sprintf('%s (%s), %s - step %d, shaded band is the %.1f deg guard sector', ...
      array.name, array.component, scenario.behaviour, step_index, guard_deg), ...
      'Color', FOREGROUND, 'FontWeight', 'normal', 'Interpreter', 'none');

xlim(axes_handle, [theta_deg(1) theta_deg(end)]);
ylim(axes_handle, y_limits);
grid(axes_handle, 'on');
set(axes_handle, 'GridColor', [0.8 0.8 0.8], 'Layer', 'top');

% ────────────────────────── SAVE ──────────────────────────────────

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

figure_path = fullfile(results_dir, sprintf('pattern_%s_%s.png', ...
    matlab.lang.makeValidName(array.name), scenario.behaviour));

% Explicit background: the session theme otherwise leaks into the export.
exportgraphics(figure_handle, figure_path, 'BackgroundColor', 'white', 'Resolution', 150);
close(figure_handle);

fprintf('saved %s\n', figure_path);
end
