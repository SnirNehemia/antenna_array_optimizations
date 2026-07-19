function output_dir = run_antijam(config_path)
% RUN_ANTIJAM  Closed-loop anti-jam campaign: scenarios x algorithms -> report.
%
%   output_dir = RUN_ANTIJAM(config_path)
%
%   Top-level driver (analogue of run_optimization for this milestone):
%     1. Load config (read_config_yaml); validate required keys of the
%        antijam / sim / adapt / agent sections — missing key -> error.
%     2. Load element patterns (load_element_patterns / stack_component per the
%        Milestone-1 polarization convention, including 'total').
%     3. Extract the 1-D azimuth cut per antijam.cut_type:
%          'phi_cut'   : E(:, nearest theta = cut_theta_deg, :), axis = phi_deg.
%          'theta_cut' : principal plane (phi = cut_phi_deg / +180 half-planes
%                        combined, complex analogue of principal_plane_cut),
%                        axis = -90..90 deg.
%     4. Build / load the codebook (agent_codebook_build, cached under
%        results/antijam/codebook_cache.mat) when the bandit is enabled.
%     5. For each scenario in antijam.scenarios x each algorithm in
%        antijam.algorithms x each Monte Carlo seed:
%          mode = 'C' for {'oracle', 'lcmv'}, 'S' for {'spsa', 'bandit'};
%          run the sim_engine_init / sim_engine_step loop with the algorithm's
%          <alg>_init / <alg>_update pair; log per-step results (run_log).
%     6. kpi_evaluate per run; aggregate; write KPI table (CSV + txt),
%        plot_antijam_report figures, and a config snapshot to a timestamped
%        results/antijam/<ts>/ folder.
%
%   Inputs:
%       config_path : path to config.yaml. Default: <repo>/config.yaml.
%
%   Outputs:
%       output_dir : the timestamped results directory that was created.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

error('run_antijam:NotImplemented', ...
    'P0 stub — assembled incrementally in P1-P6 (antijam_milestone_plan.md).');
end
