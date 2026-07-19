function [obs, sim_state] = sim_engine_step(sim_state, w)
% SIM_ENGINE_STEP  One closed-loop step: apply weights, observe SINR (+snapshots).
%
%   [obs, sim_state] = SIM_ENGINE_STEP(sim_state, w)
%
%   THE ONLY observation channel available to adaptation algorithms. Advances
%   the scenario one step (k -> k+1) and evaluates the applied weights against
%   the true jammer state, which stays hidden inside sim_state.
%
%   SINR (pattern-level, per Section 2 of antijam_milestone_plan.md):
%       P_sig = sigma_s^2 * sum_c |w' * e_c(theta_s)|^2 / n_comp
%       P_jam = sigma_j^2 * sum_c |w' * e_c(theta_j)|^2 / n_comp   (0 when off)
%       SINR  = P_sig / (P_jam + sigma_n^2 * norm(w)^2)
%   where c runs over the 1 or 2 polarization components (equal power split).
%
%   Snapshot model (Mode C):
%       x(t) = sum_c e_c(theta_s) s_c(t) + sum_c e_c(theta_j) j_c(t) + n(t)
%   with s_c, j_c, n i.i.d. circular complex Gaussian drawn from
%   sim_state.stream as sqrt(sigma^2/2) * (randn + 1i*randn).
%
%   Inputs:
%       sim_state : struct from sim_engine_init (advanced state returned).
%       w         : (N_el x 1) complex weights to apply this step.
%
%   Outputs:
%       obs : struct with fields
%           sinr_db   : scalar, output SINR of w at the current true state [dB].
%           snapshots : (N_el x K) complex, K = sim_state.n_snapshots — Mode C
%                       only; [] in Mode S. Mode-S algorithms MUST NOT touch it.
%       sim_state : updated state (k advanced; stream advances in place — it is
%                   a handle object).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P1].

k = sim_state.k + 1;
scn = sim_state.scenario;
if k > numel(scn.t_s)
    error('sim_engine_step:EndOfScenario', ...
        'Scenario ''%s'' has %d steps; step %d requested.', ...
        scn.id, numel(scn.t_s), k);
end
sim_state.k = k;

w = w(:);
n_el   = size(sim_state.E1, 1);
n_comp = 1 + double(~isempty(sim_state.E2));

% True jammer steering column(s) at the current step.
idx_j = nearest_index(sim_state.angle_deg(:), scn.theta_j_deg(k));
e_j   = sim_state.E1(:, idx_j);
if n_comp == 2
    e_j = [e_j, sim_state.E2(:, idx_j)];
end
if scn.jammer_on(k)
    sigma_j_sq = 10^(scn.jn_ratio_db(k) / 10);
else
    sigma_j_sq = 0;
end

% ── Pattern-level SINR ────────────────────────────────────────────
p_sig   = sim_state.sigma_s_sq * sum(abs(w' * sim_state.e_s).^2) / n_comp;
p_jam   = sigma_j_sq           * sum(abs(w' * e_j).^2)           / n_comp;
p_noise = sim_state.sigma_n_sq * real(w' * w);
obs = struct('sinr_db', 10 * log10(p_sig / (p_jam + p_noise)), 'snapshots', []);

% ── Mode C snapshots ──────────────────────────────────────────────
if strcmp(sim_state.mode, 'C')
    st = sim_state.stream;   % handle object: draws advance it in place
    K  = sim_state.n_snapshots;
    x  = sim_state.e_s * cgauss(st, n_comp, K, sim_state.sigma_s_sq / n_comp);
    if sigma_j_sq > 0
        x = x + e_j * cgauss(st, n_comp, K, sigma_j_sq / n_comp);
    end
    x = x + cgauss(st, n_el, K, sim_state.sigma_n_sq);
    obs.snapshots = x;
end
end


function x = cgauss(stream, m, n, sigma_sq)
% (m x n) i.i.d. circular complex Gaussian samples with variance sigma_sq.
x = sqrt(sigma_sq / 2) * (randn(stream, m, n) + 1i * randn(stream, m, n));
end
