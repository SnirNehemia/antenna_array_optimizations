function log = closed_loop_run(alg, stack1, stack2, theta_deg, phi_deg, scn, aj, sim_cfg, config, codebook)
% CLOSED_LOOP_RUN  One closed-loop run of an algorithm against a scenario.
%
%   log = CLOSED_LOOP_RUN(alg, stack1, stack2, theta_deg, phi_deg, scn, aj, ...
%                         sim_cfg, config, codebook)
%
%   Modes are fixed per algorithm: oracle, lcmv -> Mode C; spsa, bandit ->
%   Mode S. The oracle applies the perfect-knowledge LCMV weights with a
%   one-step lag (R of step k -> w applied at k+1); SPSA warm-starts from the
%   quiescent MVDR beam. Shared by run_antijam and run_jammer_demo.
%
%   [P7]: operates on the full 2-D (theta, phi) far-field grid — no 1-D cut.
%
%   Inputs:
%       alg            : 'oracle' | 'lcmv' | 'spsa' | 'bandit'.
%       stack1, stack2 : (N_el x N_theta x N_phi) far-field stacks (stack2 may
%                        be [] for single-component operation).
%       theta_deg      : (1 x N_theta) elevation grid [deg].
%       phi_deg        : (1 x N_phi) azimuth grid [deg].
%       scn            : struct from sim_scenario.
%       aj             : antijam config section.
%       sim_cfg        : sim config section (seed set per run).
%       config         : full parsed config (adapt / agent sections).
%       codebook       : struct from agent_codebook_build ([] unless alg = bandit).
%
%   Outputs:
%       log : struct with sinr_db (1 x T), W (N_el x T), arm_index (1 x T,
%             NaN except bandit), grid (E1/E2/theta_deg/phi_deg/e_s), and for
%             bandit runs null_center_deg. Caller adds oracle_sinr_db.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6, P7].

mode = 'S';
if any(strcmp(alg, {'oracle', 'lcmv'})), mode = 'C'; end
st   = sim_engine_init(stack1, stack2, theta_deg, phi_deg, scn, aj, sim_cfg, mode);
n_el = size(stack1, 1);
T    = numel(scn.t_s);

switch alg
    case 'oracle'
        w = adapt_lcmv(eye(n_el), st.e_s, 0);           % quiescent start
    case 'lcmv'
        trk = adapt_tracking_init(config.adapt, st.e_s, n_el);
        w = trk.w;
    case 'spsa'
        % Warm start from the quiescent MVDR beam — the uniform start sits
        % ~10-20 dB deeper on the real array and dominates the probe budget.
        w_q = adapt_lcmv(eye(n_el), st.e_s, 0);
        sp = adapt_spsa_init(config.adapt.spsa, w_q, sim_cfg.seed + 5000);
        w = sp.w;
    case 'bandit'
        bd = agent_bandit_init(config.agent, codebook, sim_cfg.seed + 6000);
        w = bd.w;
    otherwise
        error('closed_loop_run:BadAlgorithm', 'Unknown algorithm ''%s''.', alg);
end

log = struct();
log.sinr_db   = zeros(1, T);
log.W         = complex(zeros(n_el, T));
log.arm_index = NaN(1, T);
for k = 1:T
    log.W(:, k) = w;
    [obs, st] = sim_engine_step(st, w);
    log.sinr_db(k) = obs.sinr_db;
    switch alg
        case 'oracle'
            w = adapt_lcmv(sim_analytic_covariance(st), st.e_s, 0);
        case 'lcmv'
            [w, trk] = adapt_tracking_update(trk, obs);
        case 'spsa'
            [w, sp] = adapt_spsa_update(sp, obs);
        case 'bandit'
            log.arm_index(k) = bd.arm_index;
            [w, bd] = agent_bandit_update(bd, obs);
    end
end
log.grid = struct('E1', st.E1, 'E2', st.E2, 'theta_deg', theta_deg, ...
    'phi_deg', phi_deg, 'e_s', st.e_s);
if strcmp(alg, 'bandit')
    log.null_center_deg = codebook.null_center_deg;   % for the arm heatmap
end
end
