function R = sim_analytic_covariance(sim_state)
% SIM_ANALYTIC_COVARIANCE  True interference-plus-noise covariance (oracle only).
%
%   R = SIM_ANALYTIC_COVARIANCE(sim_state)
%
%   Analytic R at the CURRENT scenario step (the one last advanced to by
%   sim_engine_step), from ground-truth jammer state:
%       R = sigma_j^2/n_comp * e_j * e_j' + sigma_n^2 * I
%   with e_j the (N_el x n_comp) jammer steering column(s) — the matrix product
%   sums the per-component rank-1 terms. The jammer term is absent while the
%   jammer is off. NOTE: interference-plus-noise only (no desired-signal term),
%   i.e. the R of the MVDR oracle formula. Consumed by the oracle algorithm and
%   by kpi_evaluate (oracle-gap KPI). NEVER available to adapt_* / agent_*.
%
%   Inputs:
%       sim_state : struct from sim_engine_init, stepped at least once.
%
%   Outputs:
%       R : (N_el x N_el) complex Hermitian covariance.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P1].

if sim_state.k < 1
    error('sim_analytic_covariance:NotStepped', ...
        'sim_engine_step must run at least once before querying analytic R.');
end

scn    = sim_state.scenario;
k      = sim_state.k;
n_el   = size(sim_state.E1, 1);
n_comp = 1 + double(~isempty(sim_state.E2));

R = sim_state.sigma_n_sq * eye(n_el);
if scn.jammer_on(k)
    [it_j, ip_j] = nearest_index_2d(sim_state.theta_deg, sim_state.phi_deg, ...
        scn.theta_j_deg(k), scn.phi_j_deg(k));
    idx_j = (ip_j - 1) * numel(sim_state.theta_deg) + it_j;
    e_j   = sim_state.E1(:, idx_j);
    if n_comp == 2
        e_j = [e_j, sim_state.E2(:, idx_j)];
    end
    sigma_j_sq = 10^(scn.jn_ratio_db(k) / 10);
    R = R + (sigma_j_sq / n_comp) * (e_j * e_j');
end
end
