function tests = test_antijam_tracking
% TEST_ANTIJAM_TRACKING  P2 gates for the Mode C covariance tracker
% (antijam_milestone_plan.md Section 4): steady-state oracle gap < 1 dB on
% the S1-S5 suite, recovery after a 10 deg jump <= 25 updates, SINR
% availability >= 95% on the drifting scenario. Tuning per the P2 sweep:
% forgetting_lambda = 0.98, diagonal_loading_db = +10 (MPDR snapshots —
% signal present — need heavy loading against self-nulling).
%
% The SINR threshold here is 5 dB, not the config default 10 dB: the 8-element
% unit-gain toy ULA tops out at SINR_max = sigma_s^2 * N_el ~ 9 dB, so the
% operational threshold is array-dependent (customer-set for the real array).
%
% [P7]: native 2-D engine — toy ULA lives on physical theta domain [0, 180],
% boresight (theta_s_deg) at the pole (0), phi fixed at 0.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));

aj = struct('theta_s_deg', 0.0, 'phi_s_deg', 0.0, 'guard_deg', 10.0, ...
            'jn_ratio_db', 20.0, 'sigma_s_db', 0.0);
sc = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, 'seed', 1234);
acfg = struct('forgetting_lambda', 0.98, 'diagonal_loading_db', 10);
theta_deg = linspace(0, 180, 361);
phi_deg   = 0;
E    = exp(1i * pi * (0:7).' * sind(theta_deg));   % half-wavelength 8-el ULA
stack = reshape(E, 8, 361, 1);

testCase.TestData.aj = aj;
testCase.TestData.sc = sc;
testCase.TestData.acfg = acfg;
testCase.TestData.theta_deg = theta_deg;
testCase.TestData.phi_deg = phi_deg;
testCase.TestData.E = E;
testCase.TestData.stack = stack;
testCase.TestData.thr_db = 5.0;
testCase.TestData.suite = { ...
    struct('id', 'S1', 'motion', 'static', 'power', 'constant'), ...
    struct('id', 'S2', 'motion', 'drift', 'theta_drift_deg_per_s', 2.0, ...
           'phi_drift_deg_per_s', 0.0, 'power', 'constant'), ...
    struct('id', 'S3', 'motion', 'drift', 'theta_drift_deg_per_s', 2.0, ...
           'phi_drift_deg_per_s', 0.0, ...
           'theta_jump_deg', 10.0, 'jump_time_s', 30.0, 'power', 'constant'), ...
    struct('id', 'S4', 'motion', 'static', 'power', 'step', ...
           'power_step_db', 10.0, 'step_time_s', 30.0), ...
    struct('id', 'S5', 'motion', 'drift', 'theta_drift_deg_per_s', 2.0, ...
           'phi_drift_deg_per_s', 0.0, 'power', 'onoff', ...
           'duty_cycle', 0.5, 'toggle_period_s', 10.0)};
end


% ────────────────────────── HELPERS ───────────────────────────────

function log = track_loop(td, scn)
% Honest closed loop: w applied at step k comes from the update on obs(k-1).
% Oracle SINR via the MVDR max-SINR identity sigma_s^2 * e_s' R^-1 e_s.
st  = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, td.sc, 'C');
trk = adapt_tracking_init(td.acfg, st.e_s, size(td.E, 1));
T   = numel(scn.t_s);
log = struct('sinr_db', zeros(1, T), 'oracle_sinr_db', zeros(1, T));
w   = trk.w;
for k = 1:T
    [obs, st] = sim_engine_step(st, w);
    log.sinr_db(k) = obs.sinr_db;
    R = sim_analytic_covariance(st);
    log.oracle_sinr_db(k) = 10 * log10(st.sigma_s_sq * real(st.e_s' * (R \ st.e_s)));
    [w, trk] = adapt_tracking_update(trk, obs);
end
end


function ss = steady_state_mask(scn, n_skip_start, n_skip_event)
% Steady-state step indices: skip the initial convergence transient and a
% window after every scenario event (jump / power_step / turn_on).
T  = numel(scn.t_s);
ss = true(1, T);
ss(1:min(n_skip_start, T)) = false;
for i = 1:numel(scn.events)
    k_ev = find(scn.t_s >= scn.events{i}.t_s, 1);
    ss(k_ev:min(k_ev + n_skip_event, T)) = false;
end
ss = find(ss);
end


% ── P2 gate 1: steady-state oracle gap < 1 dB on the S1-S5 suite ──

function test_oracle_gap_scenario_suite(testCase)
td = testCase.TestData;
for i = 1:numel(td.suite)
    scn = sim_scenario(td.suite{i}, td.aj, td.sc);
    log = track_loop(td, scn);
    ss  = steady_state_mask(scn, 50, 25);
    gap = mean(log.oracle_sinr_db(ss) - log.sinr_db(ss));
    verifyLessThan(testCase, gap, 1.0, ...
        sprintf('%s: steady-state oracle gap %.3f dB', scn.id, gap));
end
end


% ── P2 gate 2: recovery after a 10 deg jump <= 25 updates ─────────

function test_recovery_after_jump(testCase)
td  = testCase.TestData;
scn = sim_scenario(td.suite{3}, td.aj, td.sc);          % S3
log = track_loop(td, scn);
k_jump = find(scn.t_s >= 30.0, 1);
% Above threshold going into the jump, back above it within 25 updates.
verifyGreaterThanOrEqual(testCase, log.sinr_db(k_jump - 1), td.thr_db);
rec = find(log.sinr_db(k_jump:end) >= td.thr_db, 1) - 1;
verifyNotEmpty(testCase, rec, 'never recovered after the jump');
verifyLessThanOrEqual(testCase, rec, 25, ...
    sprintf('recovery took %d covariance updates', rec));
end


% ── P2 gate 3: availability >= 95% on the drifting scenario ───────

function test_availability_drifting(testCase)
td  = testCase.TestData;
scn = sim_scenario(td.suite{2}, td.aj, td.sc);          % S2
log = track_loop(td, scn);
avail = mean(log.sinr_db >= td.thr_db);                 % full timeline, no grace
verifyGreaterThanOrEqual(testCase, avail, 0.95, ...
    sprintf('SINR availability %.4f', avail));
end


% ── Contract: tracker rejects Mode S observations ─────────────────

function test_tracker_requires_mode_c(testCase)
td   = testCase.TestData;
trk  = adapt_tracking_init(td.acfg, td.E(:, 1), 8);
obs_s = struct('sinr_db', 0.0, 'snapshots', []);
verifyError(testCase, @() adapt_tracking_update(trk, obs_s), ...
    'adapt_tracking_update:NoSnapshots');
% Missing config key raises (no silent defaults).
verifyError(testCase, ...
    @() adapt_tracking_init(struct('forgetting_lambda', 0.98), td.E(:, 1), 8), ...
    'adapt_tracking_init:MissingKey');
end
