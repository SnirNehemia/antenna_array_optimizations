function plot_comparison(all_scenario_results, theta_deg, n_side, cfg, output_path)
% PLOT_COMPARISON  Classical-vs-optimised comparison figure (headless).
%
%   Mirrors scripts/compare_classical.py:plot_comparison. One row per scenario:
%   left = absolute-gain (dBi) principal-plane pattern overlay for every
%   technique; right = per-directive min/mean/max whiskers across techniques.
%
%   Inputs:
%       all_scenario_results : cell array; each element a struct with fields
%                              .scenario (struct) and .results (struct keyed by
%                              technique name, as returned by run_scenario).
%       theta_deg            : one-sided elevation grid [deg].
%       n_side               : array dimension (figure title).
%       cfg                  : config struct (d_over_lambda, element_source,
%                              plot_dynamic_range_db).
%       output_path          : PNG file path.
%
%   Part of: Antenna Array Pattern Optimization Tool.

CLASSICAL = {'uniform', 'hamming', 'hanning', 'kaiser_b3', 'kaiser_b6', ...
             'chebyshev_25dB', 'chebyshev_40dB', 'taylor_25dB'};
TECHS = [CLASSICAL, {'optimized'}];

theta_display = principal_plane_theta_axis(theta_deg);
d_over_lambda = cfg_get(cfg, 'd_over_lambda', 0.5);
dyn           = cfg_get(cfg, 'plot_dynamic_range_db', 60);
if strcmp(cfg_get(cfg, 'element_source', 'synthetic'), 'folder')
    elem_label = 'CST elements';
else
    elem_label = 'isotropic elements';
end

n_rows = numel(all_scenario_results);
fig = figure('Visible', 'off', 'Position', [50 50 1500 450 * max(n_rows, 1)]);
tl  = tiledlayout(fig, n_rows, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Classical tapering/steering vs optimiser - %dx%d URA (d=%.2f \\lambda, %s)', ...
    n_side, n_side, d_over_lambda, elem_label));

for r = 1:n_rows
    scenario = all_scenario_results{r}.scenario;
    results  = all_scenario_results{r}.results;

    % ── Left: pattern overlay (dBi) ──────────────────────────────
    ax_pat = nexttile(tl);
    hold(ax_pat, 'on');
    max_peak_dbi = -inf;
    for ti = 1:numel(TECHS)
        nm = TECHS{ti};
        if ~isfield(results, nm), continue; end
        cut_vm  = results.(nm).cut_vm;
        cut_dbi = 20 * log10(max(cut_vm(:), 1e-30));
        max_peak_dbi = max(max_peak_dbi, max(cut_dbi));
        if strcmp(nm, 'optimized')
            plot(ax_pat, theta_display, cut_dbi, 'm-', 'LineWidth', 2.2, 'DisplayName', nm);
        else
            plot(ax_pat, theta_display, cut_dbi, 'LineWidth', 1.0, 'DisplayName', nm);
        end
    end
    floor_dbi = max_peak_dbi - dyn;
    ylim(ax_pat, [floor_dbi, max_peak_dbi + 2]);
    xlim(ax_pat, [theta_display(1), theta_display(end)]);
    xlabel(ax_pat, 'theta display (deg)'); ylabel(ax_pat, 'Gain (dBi)');
    title(ax_pat, sprintf('%s - pattern cut', scenario_label(scenario)));
    grid(ax_pat, 'on');
    legend(ax_pat, 'Location', 'southoutside', 'NumColumns', 3, 'FontSize', 6);

    % ── Right: per-directive whiskers ────────────────────────────
    ax_w = nexttile(tl);
    hold(ax_w, 'on');
    directives = scenario.directives;
    n_dir = numel(directives);
    present = TECHS(cellfun(@(nm) isfield(results, nm), TECHS));
    x = 1:numel(present);
    offsets = linspace(-0.18, 0.18, max(n_dir, 1));
    for di = 1:n_dir
        for xi = 1:numel(present)
            dm   = results.(present{xi}).directive_metrics(di);
            peak = dm.peak_dbi;
            lo   = dm.min_db  + peak;
            hi   = dm.max_db  + peak;
            mu   = dm.mean_db + peak;
            xc   = x(xi) + offsets(di);
            if strcmp(present{xi}, 'optimized'), col = 'm'; else, col = [0.27 0.51 0.71]; end
            plot(ax_w, [xc xc], [lo hi], '-', 'Color', col, 'LineWidth', 1.4);
            plot(ax_w, xc, mu, 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 4);
        end
    end
    set(ax_w, 'XTick', x, 'XTickLabel', present, 'XTickLabelRotation', 35, 'FontSize', 7);
    ylabel(ax_w, 'Gain (dBi)'); grid(ax_w, 'on');
    title(ax_w, 'Directive windows: [min,max] + mean dot');
end

exportgraphics(fig, output_path, 'Resolution', 150);
close(fig);
end


function s = scenario_label(scenario)
if isfield(scenario, 'label') && ~isempty(scenario.label)
    s = scenario.label;
else
    s = 'scenario';
end
end


function val = cfg_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end
