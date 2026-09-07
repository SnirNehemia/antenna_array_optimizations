function m = kpi_sweep_metrics(run_log, scenario, antijam_config, stack1, stack2, ...
                               theta_deg, phi_deg, ss_start_frac, full_kpi, ...
                               dir_ref_dbi, track_tol_db)
% KPI_SWEEP_METRICS  Scalar performance metrics for one closed-loop run.
%
%   m = KPI_SWEEP_METRICS(run_log, scenario, antijam_config, stack1, stack2, ...
%                         theta_deg, phi_deg, ss_start_frac, full_kpi, ...
%                         dir_ref_dbi, track_tol_db)
%
%   [P12] Extracted verbatim from run_amplitude_sweep_script's local
%   sweep_metrics so the amplitude sweep and the acceptance grid score runs with
%   ONE definition — two drivers each carrying a copy of a metric is how a
%   verdict disagreement like P9/P11's becomes possible. Only addition: the
%   track_score_pct headline (below). Every pre-existing field is unchanged,
%   including its dB-domain averaging.
%
%   Definitions mirror the authoritative ones verbatim:
%       availability / dead time  -> plot_mode_c_comparison
%       oracle gap / recovery     -> kpi_evaluate
%       directivity toward target -> compute_directivity_trace
%   They are recomputed here rather than taken from kpi_evaluate because that
%   function's null-pointing KPI scans the whole far-field grid for local minima
%   at every step (~95% of a cell's cost) and the sweeps do not map it. Pass
%   full_kpi = true to route through kpi_evaluate instead.
%
%   Inputs:
%       run_log        : struct from closed_loop_run, with oracle_sinr_db added
%                        by the caller. Required: sinr_db, oracle_sinr_db, W, grid.
%       scenario       : struct from sim_scenario (needs t_s, events).
%       antijam_config : struct. Required: sinr_min_db.
%       stack1, stack2 : (N_el x N_theta x N_phi) far-field stacks; stack2 is []
%                        for single-component operation. Units: V/m.
%       theta_deg, phi_deg : grid axes. Units: degrees.
%       ss_start_frac  : steady-state window start as a fraction of run duration
%                        (0.5 = last half). Units: dimensionless.
%       full_kpi       : logical; route through kpi_evaluate (slow) when true.
%       dir_ref_dbi    : quiescent-beam directivity toward the target FOR THIS
%                        ARRAY, from kpi_quiescent_directivity. Units: dBi.
%       track_tol_db   : oracle-tracking tolerance. Units: dB.
%
%   Outputs:
%       m : struct of scalar metrics; field list is kpi_sweep_metric_names().
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P11, P12].

thr  = antijam_config.sinr_min_db;
sinr = run_log.sinr_db;
t    = scenario.t_s;
dt   = t(2) - t(1);
ss   = t >= ss_start_frac * t(end);            % steady-state mask

if full_kpi
    k = kpi_evaluate(run_log, scenario, antijam_config);
    m.availability_pct    = 100 * k.availability;
    m.oracle_gap_mean_db  = k.oracle_gap_mean_db;
    m.recovery_mean_steps = mean(k.recovery_steps, 'omitnan');
    m.null_err_median_deg = median(k.null_pointing_err_deg( ...
        ~isnan(k.null_pointing_err_deg)));
else
    m.availability_pct    = 100 * mean(sinr >= thr);
    m.oracle_gap_mean_db  = mean(run_log.oracle_sinr_db - sinr);
    m.recovery_mean_steps = recovery_mean(sinr, scenario, thr);
end
m.dead_time_s      = dt * sum(sinr < thr);
% NOTE: means are taken in dB (a geometric mean in linear power), matching
% kpi_evaluate's oracle_gap_mean_db convention.
m.sinr_mean_db     = mean(sinr);
m.sinr_ss_db       = mean(sinr(ss));
m.oracle_gap_ss_db = mean(run_log.oracle_sinr_db(ss) - sinr(ss));

% ── [P12] Headline metric: oracle-tracking score ──────────────────
% Fraction of the run spent within track_tol_db of the perfect-knowledge LCMV.
% Normalized against the best ACHIEVABLE rather than a fixed threshold, which is
% what makes cells at different amplitudes, angles and arrays directly
% comparable: availability cannot rank algorithms across easy cells because it
% saturates at 100% on all of them. Slow convergence, a poor steady-state gap
% and slow recovery all cost score. The oracle scores exactly 100 by
% construction, so every grid carries its own sanity anchor.
gap = run_log.oracle_sinr_db - sinr;
m.track_score_pct    = 100 * mean(gap <= track_tol_db);
m.track_score_ss_pct = 100 * mean(gap(ss) <= track_tol_db);

dir_s = compute_directivity_trace(run_log, stack1, stack2, theta_deg, phi_deg);
m.dir_s_dbi_mean = mean(dir_s);
m.dir_s_dbi_ss   = mean(dir_s(ss));
% Directivity LOSS against the quiescent beam. Negative = the adapted beam
% points less energy at the target than doing nothing would have. Large negative
% values with a healthy SINR are desired-signal cancellation. [P11] measured
% this at r = 0.94-0.996 against oracle_gap_ss_db, so it is retained as a
% diagnostic rather than promoted to a headline.
m.dir_loss_db_ss = m.dir_s_dbi_ss - dir_ref_dbi;
end


function r = recovery_mean(sinr, scenario, thr)
% Mean steps from each scenario event until SINR re-crosses the threshold
% (kpi_evaluate section 2). NaN when the scenario defines no events.
n_ev = numel(scenario.events);
if n_ev == 0
    r = NaN;
    return
end
rec = NaN(1, n_ev);
for i = 1:n_ev
    k0 = find(scenario.t_s >= scenario.events{i}.t_s, 1);
    v  = find(sinr(k0:end) >= thr, 1) - 1;
    if ~isempty(v), rec(i) = v; end
end
r = mean(rec, 'omitnan');
end
