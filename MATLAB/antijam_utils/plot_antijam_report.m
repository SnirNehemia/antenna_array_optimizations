function plot_antijam_report(runs, config, output_dir)
% PLOT_ANTIJAM_REPORT  Headline figures for the anti-jam evaluation campaign.
%
%   PLOT_ANTIJAM_REPORT(runs, config, output_dir)
%
%   Saves to output_dir (one PNG per figure):
%       sinr_<scn>.png       SINR timelines (first seed), all algorithms +
%                            oracle + threshold band, jammer-on shading.
%       recovery_hist.png    Recovery-time histograms per algorithm
%                            (all events, all scenarios & seeds).
%       arm_heatmap_<scn>.png  Bandit arm-center track vs true theta_j(t).
%       oracle_gap_bars.png  Mean oracle gap per scenario x algorithm.
%       patterns_<scn>.png   Cut-pattern snapshots around scenario events.
%       regret.png           Cumulative oracle-gap (regret) curves, Mode S.
%       null_lifecycle_<scn>.png  For 'window'-power scenarios (S6): SINR,
%                            peak-gain penalty, null depth at theta_j, and
%                            bandit arm index through on -> off transitions.
%
%   Inputs:
%       runs       : cell array of run records from run_antijam (fields
%                    scenario_id, algorithm, seed, run_log, kpi, scenario).
%       config     : full parsed config (antijam / sim / adapt / agent).
%       output_dir : destination folder (exists).
%
%   Outputs: none (figure files on disk).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

aj   = config.antijam;
thr  = aj.sinr_min_db;
scn_ids = unique(cellfun(@(r) {r.scenario_id}, runs), 'stable');
algs    = unique(cellfun(@(r) {r.algorithm}, runs), 'stable');
cols    = lines(numel(algs));
seed0   = min(cellfun(@(r) r.seed, runs));

% ── 1. SINR timelines per scenario (first seed) ───────────────────
for s = 1:numel(scn_ids)
    f = figure('Visible', 'off', 'Position', [50 50 950 420]);
    hold on; grid on;
    scn = pick(runs, scn_ids{s}, algs{1}, seed0).scenario;
    shade_jammer_on(scn);
    for a = 1:numel(algs)
        r = pick(runs, scn_ids{s}, algs{a}, seed0);
        plot(scn.t_s, movmean(r.run_log.sinr_db, 11), 'Color', cols(a, :), ...
            'LineWidth', 1.1, 'DisplayName', algs{a});
    end
    yline(thr, 'r--', sprintf('threshold %.0f dB', thr), 'DisplayName', 'threshold');
    xlabel('t [s]'); ylabel('SINR [dB] (11-step moving mean)');
    title(sprintf('%s — SINR timeline (seed %d)', scn_ids{s}, seed0));
    legend('Location', 'southeast');
    save_png(f, output_dir, sprintf('sinr_%s.png', scn_ids{s}));
end

% ── 2. Recovery histograms per algorithm ──────────────────────────
f = figure('Visible', 'off', 'Position', [50 50 950 300]);
for a = 1:numel(algs)
    recs = [];
    for i = 1:numel(runs)
        if strcmp(runs{i}.algorithm, algs{a})
            recs = [recs, runs{i}.kpi.recovery_steps]; %#ok<AGROW>
        end
    end
    recs = recs(~isnan(recs));
    subplot(1, numel(algs), a);
    histogram(recs, 'BinWidth', 5); grid on;
    title(sprintf('%s (med %.0f)', algs{a}, median(recs)));
    xlabel('recovery [steps]');
end
save_png(f, output_dir, 'recovery_hist.png');

% ── 3. Bandit arm heatmap per scenario ────────────────────────────
if any(strcmp(algs, 'bandit'))
    for s = 1:numel(scn_ids)
        r = pick(runs, scn_ids{s}, 'bandit', seed0);
        if isempty(r), continue; end
        cb_centers = arm_centers(r);
        f = figure('Visible', 'off', 'Position', [50 50 950 380]);
        hold on; grid on;
        shade_jammer_on(r.scenario);
        played = r.run_log.arm_index;
        ok = ~isnan(played) & ~isnan(cb_centers(max(played, 1)));
        scatter(r.scenario.t_s(ok), cb_centers(played(ok)), 6, 'b', 'filled', ...
            'DisplayName', 'played arm center');
        plot(r.scenario.t_s, r.scenario.theta_j_deg, 'r-', 'LineWidth', 1.2, ...
            'DisplayName', '\theta_j (truth)');
        xlabel('t [s]'); ylabel('cut angle [deg]');
        title(sprintf('%s — bandit arm selection vs true jammer angle (seed %d)', ...
            scn_ids{s}, seed0));
        legend('Location', 'best');
        save_png(f, output_dir, sprintf('arm_heatmap_%s.png', scn_ids{s}));
    end
end

% ── 4. Oracle-gap bars ────────────────────────────────────────────
gap = zeros(numel(scn_ids), numel(algs));
for s = 1:numel(scn_ids)
    for a = 1:numel(algs)
        vals = [];
        for i = 1:numel(runs)
            if strcmp(runs{i}.scenario_id, scn_ids{s}) && strcmp(runs{i}.algorithm, algs{a})
                vals(end + 1) = runs{i}.kpi.oracle_gap_mean_db; %#ok<AGROW>
            end
        end
        gap(s, a) = mean(vals);
    end
end
f = figure('Visible', 'off', 'Position', [50 50 950 380]);
bar(gap); grid on;
set(gca, 'XTickLabel', scn_ids);
ylabel('mean oracle gap [dB]'); legend(algs, 'Location', 'northwest');
title('Oracle gap per scenario \times algorithm (Monte Carlo mean)');
save_png(f, output_dir, 'oracle_gap_bars.png');

% ── 5. Pattern snapshots around events ────────────────────────────
for s = 1:numel(scn_ids)
    r0 = pick(runs, scn_ids{s}, algs{1}, seed0);
    if isempty(r0.scenario.events), continue; end
    t_ev = r0.scenario.events{1}.t_s;
    T_end = r0.scenario.t_s(end);
    snap_t = unique(max(0, min([t_ev - 5, t_ev + 5, T_end - 1], T_end)));
    f = figure('Visible', 'off', 'Position', [50 50 1200 320]);
    plot_algs = intersect({'lcmv', 'bandit'}, algs, 'stable');
    for si = 1:numel(snap_t)
        subplot(1, numel(snap_t), si); hold on; grid on;
        k = find(r0.scenario.t_s >= snap_t(si), 1);
        for a = 1:numel(plot_algs)
            r = pick(runs, scn_ids{s}, plot_algs{a}, seed0);
            p = cut_pattern_db(r.run_log.W(:, k), r.run_log.cut);
            plot(r.run_log.cut.angle_deg, p, 'DisplayName', plot_algs{a});
        end
        xline(aj.theta_s_deg, 'g:', '\theta_s');
        if r0.scenario.jammer_on(k)
            xline(r0.scenario.theta_j_deg(k), 'r--', '\theta_j');
        end
        ylim([-70, 5]); xlabel('cut angle [deg]'); ylabel('pattern [dB]');
        title(sprintf('%s @ t = %.0f s', scn_ids{s}, snap_t(si)));
        if si == 1, legend('Location', 'southwest'); end
    end
    save_png(f, output_dir, sprintf('patterns_%s.png', scn_ids{s}));
end

% ── 6. Regret curves (Mode S algorithms) ──────────────────────────
f = figure('Visible', 'off', 'Position', [50 50 950 380]);
hold on; grid on;
mode_s = intersect({'spsa', 'bandit'}, algs, 'stable');
styles = {'-', '--'};
for s = 1:min(2, numel(scn_ids))
    for a = 1:numel(mode_s)
        r = pick(runs, scn_ids{s}, mode_s{a}, seed0);
        regret = cumsum(max(r.run_log.oracle_sinr_db - r.run_log.sinr_db, 0)) ...
                 * (r.scenario.t_s(2) - r.scenario.t_s(1));
        plot(r.scenario.t_s, regret, styles{s}, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s / %s', mode_s{a}, scn_ids{s}));
    end
end
xlabel('t [s]'); ylabel('cumulative regret [dB\cdots]');
title('Regret vs oracle (Mode S algorithms)'); legend('Location', 'northwest');
save_png(f, output_dir, 'regret.png');

% ── 7. Null lifecycle for 'window' scenarios (S6) ─────────────────
for s = 1:numel(scn_ids)
    r0 = pick(runs, scn_ids{s}, algs{1}, seed0);
    if ~is_window_scenario(r0.scenario), continue; end
    scn = r0.scenario;
    plot_algs = intersect({'lcmv', 'bandit'}, algs, 'stable');
    f = figure('Visible', 'off', 'Position', [50 50 1000 760]);
    subplot(3, 1, 1); hold on; grid on;
    shade_jammer_on(scn);
    for a = 1:numel(plot_algs)
        r = pick(runs, scn_ids{s}, plot_algs{a}, seed0);
        plot(scn.t_s, movmean(r.run_log.sinr_db, 11), 'LineWidth', 1.1, ...
            'DisplayName', plot_algs{a});
    end
    plot(scn.t_s, pick(runs, scn_ids{s}, 'oracle', seed0).run_log.sinr_db, ...
        'k:', 'DisplayName', 'oracle');
    yline(thr, 'r--', 'HandleVisibility', 'off');
    ylabel('SINR [dB]'); legend('Location', 'southeast');
    title(sprintf('%s — null lifecycle: silent \\rightarrow jammer on \\rightarrow silent (seed %d)', ...
        scn_ids{s}, seed0));

    subplot(3, 1, 2); hold on; grid on;
    shade_jammer_on(scn);
    for a = 1:numel(plot_algs)
        r = pick(runs, scn_ids{s}, plot_algs{a}, seed0);
        plot(scn.t_s, r.kpi.peak_gain_penalty_db, 'LineWidth', 1.1, ...
            'DisplayName', [plot_algs{a} ' gain penalty']);
        % Clamp for display: projection nulls are numerically exact (~-320 dB)
        % and would flatten every other trace on the axis.
        depth = max(null_depth_at_thetaj(r), -60);
        plot(scn.t_s, depth, '--', 'LineWidth', 0.9, ...
            'DisplayName', [plot_algs{a} ' depth @ \theta_j (clamped -60)']);
    end
    ylabel('[dB]'); legend('Location', 'southwest');

    subplot(3, 1, 3); hold on; grid on;
    shade_jammer_on(scn);
    if any(strcmp(algs, 'bandit'))
        r = pick(runs, scn_ids{s}, 'bandit', seed0);
        stairs(scn.t_s, r.run_log.arm_index, 'b-', 'DisplayName', 'bandit arm #');
        yline(1, 'g:', 'arm 1 = unconstrained peak', 'HandleVisibility', 'off');
    end
    xlabel('t [s]'); ylabel('arm index'); legend('Location', 'northeast');
    save_png(f, output_dir, sprintf('null_lifecycle_%s.png', scn_ids{s}));
end
end


% ────────────────────────── HELPERS ───────────────────────────────

function r = pick(runs, scn_id, alg, seed)
r = [];
for i = 1:numel(runs)
    if strcmp(runs{i}.scenario_id, scn_id) && strcmp(runs{i}.algorithm, alg) ...
            && runs{i}.seed == seed
        r = runs{i};
        return;
    end
end
end


function shade_jammer_on(scn)
% Light-red background over the jammer-on portions of the timeline.
on = scn.jammer_on(:).';
d  = diff([false, on, false]);
t0 = scn.t_s(d(1:end-1) == 1);
t1 = scn.t_s([false, d(2:end-1) == -1]);
if numel(t1) < numel(t0), t1(end + 1) = scn.t_s(end); end
yl = [-60, 60];
for i = 1:numel(t0)
    patch([t0(i), t1(i), t1(i), t0(i)], [yl(1), yl(1), yl(2), yl(2)], ...
        [1, 0.92, 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end


function p_db = cut_pattern_db(w, cut)
p = abs(w' * cut.E1).^2;
if ~isempty(cut.E2), p = p + abs(w' * cut.E2).^2; end
p_db = 10 * log10(p / max(p));
end


function centers = arm_centers(r)
% Reconstruct arm null centers from a bandit log: NaN for the unconstrained
% arm. Stored in the codebook; recovered here from the run's codebook file if
% present in the log, else approximated by played-arm indices only.
if isfield(r.run_log, 'null_center_deg')
    centers = r.run_log.null_center_deg;
else
    centers = 1:max(r.run_log.arm_index);   % fallback: index axis
end
end


function tf = is_window_scenario(scn)
% True for the off -> on -> off single-window profile (S6-style).
on = scn.jammer_on;
tf = ~on(1) && ~on(end) && nnz(diff(on)) == 2;
end


function depth = null_depth_at_thetaj(r)
scn = r.scenario;
T = numel(scn.t_s);
depth = NaN(1, T);
for k = 1:T
    if ~scn.jammer_on(k), continue; end
    p = abs(r.run_log.W(:, k)' * r.run_log.cut.E1).^2;
    if ~isempty(r.run_log.cut.E2)
        p = p + abs(r.run_log.W(:, k)' * r.run_log.cut.E2).^2;
    end
    idx = nearest_index(r.run_log.cut.angle_deg(:), scn.theta_j_deg(k));
    depth(k) = 10 * log10(p(idx) / max(p));
end
end


function save_png(f, output_dir, name)
print(f, fullfile(output_dir, name), '-dpng', '-r120');
close(f);
end
