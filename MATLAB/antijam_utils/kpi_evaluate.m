function kpi = kpi_evaluate(run_log, scenario, antijam_config)
% KPI_EVALUATE  The 5 milestone KPIs from one closed-loop run log.
%
%   kpi = KPI_EVALUATE(run_log, scenario, antijam_config)
%
%   KPI code IS allowed ground truth (scenario.theta_j_deg) — the Mode C/S
%   restriction applies to algorithms only.
%
%   Inputs:
%       run_log : struct logged by run_antijam, fields
%           t_s            : (1 x T) time grid [s].
%           sinr_db        : (1 x T) achieved SINR per step.
%           oracle_sinr_db : (1 x T) perfect-knowledge LCMV SINR per step.
%           W              : (N_el x T) applied weights per step.
%           arm_index      : (1 x T) selected arm (bandit runs; NaN otherwise).
%           pattern_db     : (N_angles x T) achieved cut pattern [dB] (for
%                            null-pointing; may be subsampled in time).
%       scenario       : struct from sim_scenario (truth + events).
%       antijam_config : struct. Required: sinr_min_db, theta_s_deg.
%
%   Outputs:
%       kpi : struct with fields
%           availability        : fraction of steps with sinr_db >= sinr_min_db.
%           recovery_steps      : (1 x n_events) steps from each scenario event
%                                 until SINR re-crosses the threshold (NaN if
%                                 never recovered).
%           null_pointing_err_deg : (1 x T) |achieved null angle - theta_j(t)|.
%           peak_gain_penalty_db  : (1 x T) main-beam gain loss vs the
%                                 jammer-free Milestone-1 optimum.
%           oracle_gap_db       : (1 x T) oracle_sinr_db - sinr_db, plus a
%                                 scalar mean in oracle_gap_mean_db.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

error('kpi_evaluate:NotImplemented', ...
    'P0 stub — implementation scheduled for P6 (antijam_milestone_plan.md).');
end
