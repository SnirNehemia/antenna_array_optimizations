function names = kpi_sweep_metric_names()
% KPI_SWEEP_METRIC_NAMES  Canonical field order of a kpi_sweep_metrics struct.
%
%   names = KPI_SWEEP_METRIC_NAMES()
%
%   [P12] Single source of truth for the metric list, so a driver's map
%   allocation, per-seed accumulation, seed statistics, CSV header and CSV rows
%   cannot drift apart — adding a metric here propagates to all of them.
%   track_score_pct leads the list because it is the headline (plan Section 5,
%   KPI 6); everything after it is a diagnostic that explains a low score.
%
%   Outputs:
%       names : (1 x N) cell array of field names, in CSV column order.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P12].

names = {'track_score_pct', 'track_score_ss_pct', ...
         'availability_pct', 'dead_time_s', 'sinr_mean_db', 'sinr_ss_db', ...
         'oracle_gap_mean_db', 'oracle_gap_ss_db', 'dir_s_dbi_mean', ...
         'dir_s_dbi_ss', 'dir_loss_db_ss', 'recovery_mean_steps'};
end
