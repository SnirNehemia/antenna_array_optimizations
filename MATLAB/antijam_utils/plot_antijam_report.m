function plot_antijam_report(runs, config, output_dir)
% PLOT_ANTIJAM_REPORT  Headline figures for the anti-jam evaluation campaign.
%
%   PLOT_ANTIJAM_REPORT(runs, config, output_dir)
%
%   Saves to output_dir (one file per figure, repo plotting conventions):
%       - SINR timeline per scenario x algorithm, sinr_min_db threshold band.
%       - Recovery-time histograms per algorithm (Monte Carlo aggregate).
%       - Arm-selection heatmap vs time with true theta_j(t) overlay (bandit).
%       - Oracle-gap bars per scenario x algorithm.
%       - Pattern snapshots (cut, dB) at selected instants with theta_s /
%         theta_j markers.
%       - Codebook coverage plot (P4 gate): worst-case null depth vs angle.
%
%   Inputs:
%       runs       : cell array of run records, each a struct with fields
%                    scenario_id, algorithm, seed, run_log, kpi (from
%                    kpi_evaluate), plus codebook meta on bandit runs.
%       config     : full parsed config (antijam / sim / adapt / agent).
%       output_dir : destination folder (created by run_antijam).
%
%   Outputs: none (figure files on disk).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

error('plot_antijam_report:NotImplemented', ...
    'P0 stub — implementation scheduled for P6 (antijam_milestone_plan.md).');
end
