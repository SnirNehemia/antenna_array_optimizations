function tests = test_antijam_sim
% TEST_ANTIJAM_SIM  P1 quantitative gates for the anti-jam simulation harness
% (antijam_milestone_plan.md Section 4): oracle null depth, covariance
% convergence, scenario drift rate + determinism, hand-computed toy SINR,
% plus observation-contract checks. Analytic/toy ground truth only — no
% Python fixtures (plan Section 8). [P7]: native 2-D (theta, phi) engine —
% toy ULA patterns live on the physical theta domain [0, 180] with boresight
% at the pole (theta_s_deg = 0, guard is then a pure theta distance
% regardless of phi), phi fixed at 0 for the 1-D-equivalent gates, plus one
% genuinely 2-D toy case.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));
testCase.TestData.dir = here;
end


% ────────────────────────── HELPERS ───────────────────────────────

function [stack, theta_deg, phi_deg] = ula_patterns(n_el, n_theta)
% Ideal half-wavelength ULA of isotropic elements over theta in [0, 180]
% (boresight at the pole, theta = 0): e_n(theta) = exp(1i * pi * n * sin(theta)).
% Single-column phi axis (degenerate 2-D grid) — analytic steering ground
% truth for the 1-D-equivalent gates.
theta_deg = linspace(0, 180, n_theta);
phi_deg   = 0;
E = exp(1i * pi * (0:n_el - 1).' * sind(theta_deg));
stack = reshape(E, n_el, n_theta, 1);
end


function [aj, sc] = base_configs()
aj = struct('theta_s_deg', 0.0, 'phi_s_deg', 0.0, 'guard_deg', 10.0, ...
            'jn_ratio_db', 20.0, 'sigma_s_db', 0.0);
sc = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, ...
            'seed', 1234);
end


% ── P1 gate 1: oracle LCMV null depth >= 40 dB at J/N = 20 dB ─────

function test_oracle_null_depth(testCase)
[aj, sc] = base_configs();
[stack, theta_deg, phi_deg] = ula_patterns(8, 361);
scn   = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), aj, sc);
state = sim_engine_init(stack, [], theta_deg, phi_deg, scn, aj, sc, 'S');
[~, state] = sim_engine_step(state, ones(8, 1));

R = sim_analytic_covariance(state);
w = adapt_lcmv(R, state.e_s, 0);

pattern_pow = abs(w' * state.E1).^2;                % (1 x N_theta) power pattern
idx_j       = nearest_index(theta_deg(:), scn.theta_j_deg(1));
null_depth_db = 10 * log10(pattern_pow(idx_j) / max(pattern_pow));
verifyLessThanOrEqual(testCase, null_depth_db, -40, ...
    sprintf('oracle null depth %.1f dB at theta_j = %.1f deg', ...
    null_depth_db, scn.theta_j_deg(1)));

% Distortionless constraint holds: unit response toward theta_s.
verifyEqual(testCase, abs(w' * state.e_s), 1.0, 'AbsTol', 1e-12);
end


% ── P1 gate 2: sample covariance -> analytic, < 5% at N = 50*N_el ─

function test_covariance_convergence(testCase)
[aj, sc] = base_configs();
n_el = 8;
[stack, theta_deg, phi_deg] = ula_patterns(n_el, 361);
scn   = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), aj, sc);
state = sim_engine_init(stack, [], theta_deg, phi_deg, scn, aj, sc, 'C');

n_snapshots = 50 * n_el;
n_steps     = ceil(n_snapshots / state.n_snapshots);
X = [];
for k = 1:n_steps
    [obs, state] = sim_engine_step(state, ones(n_el, 1));
    X = [X, obs.snapshots]; %#ok<AGROW>
end
X     = X(:, 1:n_snapshots);
R_hat = (X * X') / n_snapshots;

% Snapshots include the desired signal, so compare against the FULL analytic
% covariance E[xx'] = sigma_s^2 e_s e_s' + R_interference_plus_noise.
R_full = state.sigma_s_sq * (state.e_s * state.e_s') + sim_analytic_covariance(state);
rel_err = norm(R_hat - R_full, 'fro') / norm(R_full, 'fro');
verifyLessThan(testCase, rel_err, 0.05, ...
    sprintf('relative Frobenius error %.3f at N = %d snapshots', rel_err, n_snapshots));
end


% ── P1 gate 3: drift rate within 2%; deterministic under seed ─────

function test_scenario_drift_and_determinism(testCase)
[aj, sc] = base_configs();
drift_cfg = struct('id', 'S2', 'motion', 'drift', 'theta_drift_deg_per_s', 2.0, ...
                   'phi_drift_deg_per_s', 0.0, 'power', 'constant');
scn = sim_scenario(drift_cfg, aj, sc);

% Measured drift rate (|step| is preserved by reflection except at folds).
rate = mean(abs(diff(scn.theta_j_deg))) / sc.dt_s;
verifyLessThan(testCase, abs(rate - 2.0) / 2.0, 0.02, ...
    sprintf('measured drift rate %.4f deg/s', rate));

% Trajectory never enters the guard cap or leaves the physical theta range.
sep = angular_separation_deg(scn.theta_j_deg, scn.phi_j_deg, aj.theta_s_deg, aj.phi_s_deg);
verifyGreaterThanOrEqual(testCase, min(sep), aj.guard_deg - 1e-6);
verifyGreaterThanOrEqual(testCase, min(scn.theta_j_deg), 0 - 1e-9);
verifyLessThanOrEqual(testCase, max(scn.theta_j_deg), 180 + 1e-9);

% Deterministic under a fixed seed; different seed -> different trajectory.
scn_same = sim_scenario(drift_cfg, aj, sc);
verifyEqual(testCase, scn_same.theta_j_deg, scn.theta_j_deg);
sc_other = sc;  sc_other.seed = 4321;
scn_other = sim_scenario(drift_cfg, aj, sc_other);
verifyNotEqual(testCase, scn_other.theta_j_deg(1), scn.theta_j_deg(1));
end


% ── P1 gate 4: 2-element toy SINR matches hand computation, 1e-10 ─

function test_toy_sinr_hand_computed(testCase)
[aj, sc] = base_configs();
aj.sigma_s_db = 3.0;
% Hand-built grid: e_s = [1; 1] at (theta=0, phi=0) [boresight/pole], e_j =
% [1; -1] at (theta=40, phi=0) — single-phi-column stack.
theta_deg = [0, 40];
phi_deg   = 0;
E = [1, 1; 1, -1];
stack = reshape(E, 2, 2, 1);
% Hand-built scenario (frozen contract fields) for full control of theta_j.
scn = struct('id', 'toy', 't_s', 0, 'theta_j_deg', 40, 'phi_j_deg', 0, ...
             'jammer_on', true, 'jn_ratio_db', 20.0, 'events', {{}});
state = sim_engine_init(stack, [], theta_deg, phi_deg, scn, aj, sc, 'S');

w = [2; 1i];
[obs, ~] = sim_engine_step(state, w);
% By hand: |w' e_s|^2 = |2 - 1i|^2 = 5, |w' e_j|^2 = |2 + 1i|^2 = 5,
% ||w||^2 = 5, sigma_s^2 = 10^0.3, sigma_j^2 = 100, sigma_n^2 = 1:
%   SINR = 10^0.3 * 5 / (100 * 5 + 5)
sinr_db_hand = 10 * log10(10^0.3 * 5 / (100 * 5 + 5));
verifyEqual(testCase, obs.sinr_db, sinr_db_hand, 'RelTol', 1e-10);
end


% ── P7 gate: genuine 2-D toy case (target and jammer differ in theta AND phi)

function test_toy_sinr_2d_hand_computed(testCase)
[aj, sc] = base_configs();
aj.theta_s_deg = 30.0; aj.phi_s_deg = 90.0; aj.sigma_s_db = 3.0;
theta_deg = [30, 45];
phi_deg   = [90, 180];
% E(:, it, ip): put e_s = [1;1] at (30,90) i.e. (it=1,ip=1); e_j = [1;-1] at
% (45,180) i.e. (it=2,ip=2). Other grid points unused (zero, never indexed).
E = zeros(2, 2, 2);
E(:, 1, 1) = [1; 1];
E(:, 2, 2) = [1; -1];
scn = struct('id', 'toy2d', 't_s', 0, 'theta_j_deg', 45, 'phi_j_deg', 180, ...
             'jammer_on', true, 'jn_ratio_db', 20.0, 'events', {{}});
state = sim_engine_init(E, [], theta_deg, phi_deg, scn, aj, sc, 'S');

% e_s / e_j resolve to the correct 2-D grid points despite sharing no axis.
verifyEqual(testCase, state.e_s, [1; 1]);
w = [2; 1i];
[obs, ~] = sim_engine_step(state, w);
sinr_db_hand = 10 * log10(10^0.3 * 5 / (100 * 5 + 5));
verifyEqual(testCase, obs.sinr_db, sinr_db_hand, 'RelTol', 1e-10);

% Independent per-axis drift: theta and phi move at different rates and
% reproduce within 2% each over 60 s when far from any guard boundary.
drift_cfg = struct('id', 'D2D', 'motion', 'drift', ...
    'theta_drift_deg_per_s', 0.3, 'phi_drift_deg_per_s', -0.5, ...
    'theta_j_deg', 90.0, 'phi_j_deg', 270.0, 'power', 'constant');
scn2 = sim_scenario(drift_cfg, aj, sc);
rate_theta = mean(diff(scn2.theta_j_deg)) / sc.dt_s;
rate_phi   = mean(diff(scn2.phi_j_deg)) / sc.dt_s;   % no wrap for this small drift
verifyLessThan(testCase, abs(rate_theta - 0.3) / 0.3, 0.02);
verifyLessThan(testCase, abs(rate_phi - (-0.5)) / 0.5, 0.02);
end


% ── Fixed jammer position + per-scenario JNR override ─────────────

function test_fixed_position_and_jnr_override(testCase)
[aj, sc] = base_configs();
% Fixed position: deterministic regardless of seed, exact value.
cfg = struct('id', 'F1', 'motion', 'static', 'power', 'constant', ...
             'theta_j_deg', 40.0, 'phi_j_deg', 0.0, 'jn_ratio_db', 27.0);
scn = sim_scenario(cfg, aj, sc);
verifyEqual(testCase, scn.theta_j_deg(1), 40.0, 'AbsTol', 1e-12);
sc2 = sc;  sc2.seed = 999;
scn2 = sim_scenario(cfg, aj, sc2);
verifyEqual(testCase, scn2.theta_j_deg, scn.theta_j_deg);
% Per-scenario JNR override wins over the antijam default (20 dB).
verifyEqual(testCase, scn.jn_ratio_db(1), 27.0);
% Drift from a fixed position starts exactly there.
cfg_d = struct('id', 'F2', 'motion', 'drift', 'theta_drift_deg_per_s', 2.0, ...
               'phi_drift_deg_per_s', 0.0, 'power', 'constant', ...
               'theta_j_deg', 60.0, 'phi_j_deg', 0.0);
scn_d = sim_scenario(cfg_d, aj, sc);
verifyEqual(testCase, scn_d.theta_j_deg(1), 60.0, 'AbsTol', 1e-12);
% A position inside the guard cap raises.
cfg_bad = struct('id', 'F3', 'motion', 'static', 'power', 'constant', ...
                 'theta_j_deg', 3.0, 'phi_j_deg', 0.0);
verifyError(testCase, @() sim_scenario(cfg_bad, aj, sc), ...
    'sim_scenario:AngleInGuard');
end


% ── Observation contract: Mode C vs Mode S, error paths ───────────

function test_observation_contract(testCase)
[aj, sc] = base_configs();
[stack, theta_deg, phi_deg] = ula_patterns(4, 181);
scn = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), aj, sc);

state_c = sim_engine_init(stack, [], theta_deg, phi_deg, scn, aj, sc, 'C');
[obs_c, state_c] = sim_engine_step(state_c, ones(4, 1));
verifySize(testCase, obs_c.snapshots, [4, sc.snapshots_per_step]);

state_s = sim_engine_init(stack, [], theta_deg, phi_deg, scn, aj, sc, 'S');
[obs_s, ~] = sim_engine_step(state_s, ones(4, 1));
verifyEmpty(testCase, obs_s.snapshots);

% Same seed => identical SINR observation across modes (scalar channel is
% deterministic given the scenario; snapshots draw from the private stream).
verifyEqual(testCase, obs_s.sinr_db, obs_c.sinr_db, 'AbsTol', 1e-12);

% Error paths: covariance before stepping; stepping past scenario end.
state_fresh = sim_engine_init(stack, [], theta_deg, phi_deg, scn, aj, sc, 'S');
verifyError(testCase, @() sim_analytic_covariance(state_fresh), ...
    'sim_analytic_covariance:NotStepped');
short_scn = scn;
verifyError(testCase, @() step_past_end(stack, theta_deg, phi_deg, short_scn, aj, sc), ...
    'sim_engine_step:EndOfScenario');

% Missing required config key raises a descriptive error (no silent defaults).
aj_broken = rmfield(aj, 'guard_deg');
verifyError(testCase, ...
    @() sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), aj_broken, sc), ...
    'sim_scenario:MissingKey');
end


function step_past_end(stack, theta_deg, phi_deg, scn, aj, sc)
% Drive the engine one step beyond the scenario end.
scn.t_s = scn.t_s(1);  scn.theta_j_deg = scn.theta_j_deg(1);
scn.phi_j_deg = scn.phi_j_deg(1);
scn.jammer_on = scn.jammer_on(1);  scn.jn_ratio_db = scn.jn_ratio_db(1);
state = sim_engine_init(stack, [], theta_deg, phi_deg, scn, aj, sc, 'S');
[~, state] = sim_engine_step(state, ones(size(stack, 1), 1));
sim_engine_step(state, ones(size(stack, 1), 1));
end
