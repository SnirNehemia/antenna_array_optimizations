function tests = test_antijam_sim
% TEST_ANTIJAM_SIM  P1 quantitative gates for the anti-jam simulation harness
% (antijam_milestone_plan.md Section 4): oracle null depth, covariance
% convergence, scenario drift rate + determinism, hand-computed toy SINR,
% plus observation-contract checks. Analytic/toy ground truth only — no
% Python fixtures (plan Section 8).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));
testCase.TestData.dir = here;
end


% ────────────────────────── HELPERS ───────────────────────────────

function [E, ang] = ula_patterns(n_el, n_ang)
% Ideal half-wavelength ULA of isotropic elements over a -90..90 deg cut:
% e_n(theta) = exp(1i * pi * n * sin(theta)). Analytic steering ground truth.
ang = linspace(-90, 90, n_ang);
E   = exp(1i * pi * (0:n_el - 1).' * sind(ang));
end


function [aj, sc] = base_configs()
aj = struct('theta_s_deg', 0.0, 'guard_deg', 10.0, 'jn_ratio_db', 20.0, ...
            'sigma_s_db', 0.0, 'cut_type', 'theta_cut');
sc = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, ...
            'seed', 1234);
end


% ── P1 gate 1: oracle LCMV null depth >= 40 dB at J/N = 20 dB ─────

function test_oracle_null_depth(testCase)
[aj, sc] = base_configs();
[E, ang] = ula_patterns(8, 361);
scn   = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), aj, sc);
state = sim_engine_init(E, [], ang, scn, aj, sc, 'S');
[~, state] = sim_engine_step(state, ones(8, 1));

R = sim_analytic_covariance(state);
w = adapt_lcmv(R, state.e_s, 0);

pattern_pow = abs(w' * E).^2;                       % (1 x N_ang) power pattern
idx_j       = nearest_index(ang(:), scn.theta_j_deg(1));
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
[E, ang] = ula_patterns(n_el, 361);
scn   = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), aj, sc);
state = sim_engine_init(E, [], ang, scn, aj, sc, 'C');

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
drift_cfg = struct('id', 'S2', 'motion', 'drift', 'drift_deg_per_s', 2.0, ...
                   'power', 'constant');
scn = sim_scenario(drift_cfg, aj, sc);

% Measured drift rate (|step| is preserved by reflection except at folds).
rate = mean(abs(diff(scn.theta_j_deg))) / sc.dt_s;
verifyLessThan(testCase, abs(rate - 2.0) / 2.0, 0.02, ...
    sprintf('measured drift rate %.4f deg/s', rate));

% Trajectory never enters the guard sector or leaves the cut span.
verifyGreaterThanOrEqual(testCase, min(abs(scn.theta_j_deg - aj.theta_s_deg)), ...
    aj.guard_deg - 1e-9);
verifyLessThanOrEqual(testCase, max(abs(scn.theta_j_deg)), 90 + 1e-9);

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
% Hand-built cut: e_s = [1; 1] at 0 deg, e_j = [1; -1] at 40 deg.
E   = [1, 1; 1, -1];
ang = [0, 40];
% Hand-built scenario (frozen contract fields) for full control of theta_j.
scn = struct('id', 'toy', 't_s', 0, 'theta_j_deg', 40, 'jammer_on', true, ...
             'jn_ratio_db', 20.0, 'events', {{}});
state = sim_engine_init(E, [], ang, scn, aj, sc, 'S');

w = [2; 1i];
[obs, ~] = sim_engine_step(state, w);
% By hand: |w' e_s|^2 = |2 - 1i|^2 = 5, |w' e_j|^2 = |2 + 1i|^2 = 5,
% ||w||^2 = 5, sigma_s^2 = 10^0.3, sigma_j^2 = 100, sigma_n^2 = 1:
%   SINR = 10^0.3 * 5 / (100 * 5 + 5)
sinr_db_hand = 10 * log10(10^0.3 * 5 / (100 * 5 + 5));
verifyEqual(testCase, obs.sinr_db, sinr_db_hand, 'RelTol', 1e-10);
end


% ── Observation contract: Mode C vs Mode S, error paths ───────────

function test_observation_contract(testCase)
[aj, sc] = base_configs();
[E, ang] = ula_patterns(4, 181);
scn = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), aj, sc);

state_c = sim_engine_init(E, [], ang, scn, aj, sc, 'C');
[obs_c, state_c] = sim_engine_step(state_c, ones(4, 1));
verifySize(testCase, obs_c.snapshots, [4, sc.snapshots_per_step]);

state_s = sim_engine_init(E, [], ang, scn, aj, sc, 'S');
[obs_s, ~] = sim_engine_step(state_s, ones(4, 1));
verifyEmpty(testCase, obs_s.snapshots);

% Same seed => identical SINR observation across modes (scalar channel is
% deterministic given the scenario; snapshots draw from the private stream).
verifyEqual(testCase, obs_s.sinr_db, obs_c.sinr_db, 'AbsTol', 1e-12);

% Error paths: covariance before stepping; stepping past scenario end.
state_fresh = sim_engine_init(E, [], ang, scn, aj, sc, 'S');
verifyError(testCase, @() sim_analytic_covariance(state_fresh), ...
    'sim_analytic_covariance:NotStepped');
short_scn = scn;
verifyError(testCase, @() step_past_end(E, ang, short_scn, aj, sc), ...
    'sim_engine_step:EndOfScenario');

% Missing required config key raises a descriptive error (no silent defaults).
aj_broken = rmfield(aj, 'guard_deg');
verifyError(testCase, ...
    @() sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), aj_broken, sc), ...
    'sim_scenario:MissingKey');
end


function step_past_end(E, ang, scn, aj, sc)
% Drive the engine one step beyond the scenario end.
scn.t_s = scn.t_s(1);  scn.theta_j_deg = scn.theta_j_deg(1);
scn.jammer_on = scn.jammer_on(1);  scn.jn_ratio_db = scn.jn_ratio_db(1);
state = sim_engine_init(E, [], ang, scn, aj, sc, 'S');
[~, state] = sim_engine_step(state, ones(size(E, 1), 1));
sim_engine_step(state, ones(size(E, 1), 1));
end
