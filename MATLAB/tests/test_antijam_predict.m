function tests = test_antijam_predict
% TEST_ANTIJAM_PREDICT  P8 gates for the MUSIC + predictive Mode C nuller
% (antijam_milestone_plan.md, P8):
%   1. adapt_lcmv_null reduces to adapt_lcmv when e_null = [] (machine zero),
%      and places a deep null at the constrained direction.
%   2. adapt_music_doa recovers the jammer angle on a toy analytic covariance
%      and detects presence (on) / absence (off) via the eigengap.
%   3. On an on/off scenario the predictor is never worse than the reactive
%      lcmv baseline (availability), learns the toggle period, and estimates
%      the jammer angle accurately while it is present.
%   4. Contract: rejects Mode S observations; missing config keys raise.
%
% Toy 8-element half-wavelength ULA on the physical theta domain [0,180],
% boresight at the pole (theta_s = 0), phi fixed at 0 — same fixture family as
% test_antijam_tracking. SINR threshold 5 dB (toy SINR_max ~ 9 dB).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));

aj = struct('theta_s_deg', 0.0, 'phi_s_deg', 0.0, 'guard_deg', 10.0, ...
            'jn_ratio_db', 20.0, 'sigma_s_db', 0.0);
sc = struct('dt_s', 0.05, 'duration_s', 100.0, 'snapshots_per_step', 16, 'seed', 1234);
predict_cfg = struct('presence_gap_db', 6.0, 'buffer_len', 1024, ...
    'min_periods', 2, 'lead_steps', 6, 'doa_stride', 1);
acfg = struct('forgetting_lambda', 0.98, 'diagonal_loading_db', 10, ...
    'predict', predict_cfg);
theta_deg = linspace(0, 180, 361);                 % 0.5 deg grid
phi_deg   = 0;
E     = exp(1i * pi * (0:7).' * sind(theta_deg));  % half-wavelength 8-el ULA
stack = reshape(E, 8, 361, 1);

testCase.TestData.aj = aj;
testCase.TestData.sc = sc;
testCase.TestData.acfg = acfg;
testCase.TestData.theta_deg = theta_deg;
testCase.TestData.phi_deg = phi_deg;
testCase.TestData.E = E;
testCase.TestData.stack = stack;
testCase.TestData.thr_db = 5.0;
% On/off jammer with a 5 s off-gap (> lcmv memory ~2.5 s at lambda 0.98).
testCase.TestData.scn_onoff = struct('id', 'ONOFF', 'motion', 'static', ...
    'theta_j_deg', 30.0, 'phi_j_deg', 0.0, 'power', 'onoff', ...
    'duty_cycle', 0.5, 'toggle_period_s', 10.0);
end


% ────────────────────────── HELPERS ───────────────────────────────

function [log, prlast] = predict_loop(td, scn)
% Closed loop for the predictive nuller (w at step k from the update on obs(k-1)).
st = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, td.sc, 'C');
pr = adapt_predict_init(td.acfg, td.aj, st.e_s, td.stack, [], ...
    td.theta_deg, td.phi_deg, size(td.E, 1));
T  = numel(scn.t_s);
log = struct('sinr_db', zeros(1, T), 'doa_theta_deg', NaN(1, T), ...
    'present', false(1, T), 'period_est', NaN(1, T));
w = pr.w;
for k = 1:T
    [obs, st] = sim_engine_step(st, w);
    log.sinr_db(k) = obs.sinr_db;
    [w, pr] = adapt_predict_update(pr, obs);
    log.doa_theta_deg(k) = pr.last.theta_j_deg;
    log.present(k)       = pr.last.present;
    log.period_est(k)    = pr.last.period_est;
end
prlast = pr;
end


function log = lcmv_loop(td, scn)
% Reactive lcmv baseline (adapt_tracking) — same closed-loop convention.
st  = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, td.sc, 'C');
trk = adapt_tracking_init(td.acfg, st.e_s, size(td.E, 1));
T   = numel(scn.t_s);
log = struct('sinr_db', zeros(1, T));
w   = trk.w;
for k = 1:T
    [obs, st] = sim_engine_step(st, w);
    log.sinr_db(k) = obs.sinr_db;
    [w, trk] = adapt_tracking_update(trk, obs);
end
end


% ── Gate 1: adapt_lcmv_null correctness ───────────────────────────

function test_lcmv_null_matches_and_nulls(testCase)
td = testCase.TestData;
e_s = td.E(:, find(td.theta_deg == 0, 1));          % boresight steering
e_j = td.E(:, find(td.theta_deg == 30, 1));         % jammer steering
R   = 2 * (e_s * e_s') + 100 * (e_j * e_j') + eye(8);
load_lin = 10^(10 / 10);

% e_null = [] reproduces adapt_lcmv exactly.
w0 = adapt_lcmv_null(R, e_s, [], load_lin);
w1 = adapt_lcmv(R, e_s, load_lin);
verifyLessThan(testCase, max(abs(w0 - w1)), 1e-12, 'e_null=[] must equal adapt_lcmv');

% Explicit null: distortionless at e_s, deep null at e_j.
wn = adapt_lcmv_null(R, e_s, e_j, load_lin);
verifyEqual(testCase, abs(wn' * e_s), 1.0, 'AbsTol', 1e-9);
null_depth_db = 10 * log10(abs(wn' * e_j)^2 / abs(wn' * e_s)^2);
verifyLessThan(testCase, null_depth_db, -60, ...
    sprintf('constrained null depth %.1f dB', null_depth_db));
end


% ── Gate 2: MUSIC recovers the jammer angle + presence ────────────

function test_music_doa_recovers_angle(testCase)
td = testCase.TestData;
e_s = td.E(:, find(td.theta_deg == 0, 1));
e_j = td.E(:, find(td.theta_deg == 30, 1));
cfg = struct('theta_s_deg', 0, 'phi_s_deg', 0, 'guard_deg', 10, ...
    'presence_gap_db', 6, 'doa_stride', 1);

% Jammer ON (2 sources: signal + jammer) -> present, angle within one grid cell.
R_on = 2 * (e_s * e_s') + 100 * (e_j * e_j') + eye(8);
d_on = adapt_music_doa(R_on, td.stack, [], td.theta_deg, td.phi_deg, cfg);
verifyTrue(testCase, d_on.present, 'jammer-on not detected');
verifyLessThanOrEqual(testCase, abs(d_on.theta_j_deg - 30), 0.5, ...
    sprintf('MUSIC theta estimate %.2f (true 30)', d_on.theta_j_deg));

% Jammer OFF (signal only) -> not present.
R_off = 2 * (e_s * e_s') + eye(8);
d_off = adapt_music_doa(R_off, td.stack, [], td.theta_deg, td.phi_deg, cfg);
verifyFalse(testCase, d_off.present, 'false jammer detection when off');
end


% ── Gate 3: no regression vs lcmv + learns period + accurate DoA ──
% The performance WIN (lower dead time) is demonstrated on the real array by
% run_mode_c_demo_script; on this trivial toy ULA lcmv already reaches ~100%
% availability, so the unit gate only checks that predict does not regress and
% that its MUSIC/periodogram machinery is accurate.

function test_predict_no_regression_and_estimates(testCase)
td  = testCase.TestData;
scn = sim_scenario(td.scn_onoff, td.aj, td.sc);
[plog, ~] = predict_loop(td, scn);
llog      = lcmv_loop(td, scn);

avail_p = mean(plog.sinr_db >= td.thr_db);
avail_l = mean(llog.sinr_db >= td.thr_db);
verifyGreaterThanOrEqual(testCase, avail_p, 0.99, ...
    sprintf('predict availability %.4f too low', avail_p));
verifyGreaterThanOrEqual(testCase, avail_p, avail_l - 0.005, ...
    sprintf('predict availability %.4f regressed vs lcmv %.4f', avail_p, avail_l));

% Periodogram recovers the 10 s toggle period (200 steps) within one FFT bin.
true_period_steps = td.scn_onoff.toggle_period_s / td.sc.dt_s;
verifyLessThan(testCase, abs(plog.period_est(end) - true_period_steps) / true_period_steps, ...
    0.10, sprintf('period estimate %.1f steps (true %.0f)', ...
    plog.period_est(end), true_period_steps));

% MUSIC angle accurate while the jammer is actually transmitting (present AND
% truly on — presence lingers into the OFF gap on decaying covariance memory).
on   = plog.present & scn.jammer_on;
rmse = sqrt(mean((plog.doa_theta_deg(on) - scn.theta_j_deg(on)).^2, 'omitnan'));
verifyLessThan(testCase, rmse, 1.0, sprintf('DoA RMSE %.2f deg', rmse));
end


% ── Gate 4: Mode C contract + no silent defaults ──────────────────

function test_predict_contract(testCase)
td = testCase.TestData;
pr = adapt_predict_init(td.acfg, td.aj, td.E(:, 1), td.stack, [], ...
    td.theta_deg, td.phi_deg, 8);
obs_s = struct('sinr_db', 0.0, 'snapshots', []);
verifyError(testCase, @() adapt_predict_update(pr, obs_s), ...
    'adapt_predict_update:NoSnapshots');

% Missing adapt.predict block raises (no silent defaults).
bad = struct('forgetting_lambda', 0.98, 'diagonal_loading_db', 10);
verifyError(testCase, ...
    @() adapt_predict_init(bad, td.aj, td.E(:, 1), td.stack, [], ...
    td.theta_deg, td.phi_deg, 8), 'adapt_predict_init:MissingKey');
end


% ── Gate 5: total-pol (rank-2) beamformer reaches the true max SINR ──
% [P8 fix] A rank-2 desired signal ('total' polarization) must be received by
% the MAX-SINR beamformer (generalized Rayleigh quotient), NOT a distortionless
% constraint on every component (over-constrained, self-limiting when the two
% components differ in magnitude, e.g. co-pol vs cross-pol).

function test_total_pol_reaches_max_sinr(testCase)
td = testCase.TestData;
e1 = td.E(:, find(td.theta_deg == 30, 1));
e2 = e1 .* (1:8).';                                 % independent 2nd component
e_s = [e1, e2];
ej  = td.E(:, find(td.theta_deg == 60, 1));
R   = eye(8) + 5 * (ej * ej');                       % interference + white noise
Rs  = e1 * e1' + e2 * e2';                            % rank-2 signal covariance

% adapt_lcmv (rank-2 branch) must reach the generalized-eigenvalue max SINR.
w        = adapt_lcmv(R, e_s, 0);
sinr_ach = real(w' * Rs * w) / real(w' * R * w);
sinr_max = max(real(eig(Rs, R)));
verifyEqual(testCase, sinr_ach, sinr_max, 'RelTol', 1e-6, ...
    'rank-2 adapt_lcmv must reach the generalized-eigenvalue max SINR');

% ...and beat the old distortionless-on-both-components LCMV.
C = e_s; A = R \ C; w_old = A * ((C' * A) \ ones(2, 1));
sinr_old = real(w_old' * Rs * w_old) / real(w_old' * R * w_old);
verifyGreaterThan(testCase, sinr_ach, sinr_old, ...
    'max-SINR must beat the over-constrained 2-constraint LCMV');

% adapt_lcmv_null with e_null = [] agrees with adapt_lcmv on the rank-2 path.
w_n    = adapt_lcmv_null(R, e_s, [], 0);
sinr_n = real(w_n' * Rs * w_n) / real(w_n' * R * w_n);
verifyEqual(testCase, sinr_n, sinr_max, 'RelTol', 1e-6, ...
    'adapt_lcmv_null (e_null=[]) must match adapt_lcmv for rank-2');

% The hard null must still be exact on the rank-2 path.
w_hn = adapt_lcmv_null(R, e_s, ej, 0);
verifyLessThan(testCase, abs(w_hn' * ej) / norm(w_hn), 1e-9, ...
    'rank-2 adapt_lcmv_null must place an exact hard null');
end
