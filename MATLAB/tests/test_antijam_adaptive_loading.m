function tests = test_antijam_adaptive_loading
% TEST_ANTIJAM_ADAPTIVE_LOADING  P9 gates for the data-driven diagonal loading
% (antijam_milestone_plan.md Section 4, P9): opt-in, no regression at the
% original P2-tuned operating point, and — the actual point of the phase —
% no per-scenario re-tune needed across a wide sigma_s_db sweep (unlike the
% fixed diagonal_loading_db, whose gap blows up as sigma_s_db moves away from
% the point it was swept against).
%
% Same toy 8-element half-wavelength ULA as test_antijam_tracking.m
% (unit-gain elements, K=16 snapshots/step, lambda=0.90). Sweeping sigma_s_db
% here (jn_ratio_db fixed at 20 dB) reproduces the same regime shift that
% broke the real-array mode_c_demo (desired signal moving from far below to
% far above the jammer's eigenvalue), on a fast, reproducible fixture.
%
% [P9, 2026-08-02] Empirical calibration of the thresholds below (see
% docs/notes.md [P9] entry): on this toy fixture the adaptive scheme is
% essentially TIED with the fixed diagonal_loading_db=10 at/near the original
% calibrated point (sigma_s_db <= 10 dB; regression is a few hundredths of a
% dB, within run-to-run noise) and pulls ahead as sigma_s_db grows past it
% (measured improvement +1.16 dB at sigma_s_db=30, +3.1 dB at sigma_s_db=40) —
% smaller than the real embedded-array demo's improvement (13.6->7.9 dB gap)
% because this toy ULA's unit-gain elements and K=16 (vs the demo's K=128)
% both cap how much loading tuning alone can buy back (see the batching-fix
% note in P8's plan section: K is the dominant lever on the estimation-noise
% floor). The gate therefore checks (1) no meaningful regression anywhere in
% the sweep and (2) a real, conservatively-thresholded improvement at the
% high end — not a uniform <1 dB bound, which would not be an honest claim
% for this regime (see the P9 plan section's brute-force sweep: ~8 dB is the
% realistic ceiling, not ~1 dB, once the signal dominates the jammer).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));

aj = struct('theta_s_deg', 0.0, 'phi_s_deg', 0.0, 'guard_deg', 10.0, ...
            'jn_ratio_db', 20.0, 'sigma_s_db', 0.0);
sc = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, 'seed', 1234);
theta_deg = linspace(0, 180, 361);
phi_deg   = 0;
E    = exp(1i * pi * (0:7).' * sind(theta_deg));   % half-wavelength 8-el ULA
stack = reshape(E, 8, 361, 1);

testCase.TestData.aj = aj;
testCase.TestData.sc = sc;
testCase.TestData.theta_deg = theta_deg;
testCase.TestData.phi_deg = phi_deg;
testCase.TestData.E = E;
testCase.TestData.stack = stack;
testCase.TestData.acfg_fixed = struct('forgetting_lambda', 0.90, 'diagonal_loading_db', 10);
testCase.TestData.acfg_adapt = struct('forgetting_lambda', 0.90, ...
    'diagonal_loading_db', 10, 'loading_factor_db', 0);
testCase.TestData.scn_cfg_s1 = struct('id', 'S1', 'motion', 'static', 'power', 'constant');
end


% ────────────────────────── HELPERS ───────────────────────────────

function log = track_loop(td, aj, scn, acfg)
% Same honest closed loop as test_antijam_tracking.m's track_loop, but takes
% aj/acfg as explicit parameters so fixed vs. adaptive loading, and swept
% sigma_s_db values, can be run side by side without touching td.aj.
st  = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, aj, td.sc, 'C');
trk = adapt_tracking_init(acfg, st.e_s, size(td.E, 1));
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


function gap = steady_state_gap(td, aj, scn, acfg)
log = track_loop(td, aj, scn, acfg);
ss  = 51:numel(scn.t_s);                      % skip the initial transient
gap = mean(log.oracle_sinr_db(ss) - log.sinr_db(ss));
end


% ── Contract: absent loading_factor_db is bit-identical to pre-P9 ─

function test_adaptive_loading_off_by_default(testCase)
td = testCase.TestData;
e_s = td.E(:, 1);
trk_old = adapt_tracking_init(td.acfg_fixed, e_s, 8);
verifyFalse(testCase, trk_old.adaptive_loading);

acfg_absent = td.acfg_fixed;                  % no loading_factor_db field at all
trk_absent  = adapt_tracking_init(acfg_absent, e_s, 8);
verifyFalse(testCase, trk_absent.adaptive_loading);
verifyEqual(testCase, trk_absent.loading, trk_old.loading);

acfg_empty = td.acfg_fixed;
acfg_empty.loading_factor_db = [];             % present but empty -> also off
trk_empty = adapt_tracking_init(acfg_empty, e_s, 8);
verifyFalse(testCase, trk_empty.adaptive_loading);
end


% ── Contract: too few elements for the 2*n_comp signal+jammer split ─

function test_adaptive_loading_too_few_elements(testCase)
td = testCase.TestData;
e_s = td.E(1:2, 1);                            % N_el = 2 = n_sig -> must error
verifyError(testCase, ...
    @() adapt_tracking_init(td.acfg_adapt, e_s, 2), ...
    'adapt_tracking_init:TooFewElements');
end


% ── No regression at the original P2-tuned operating point ────────

function test_adaptive_loading_no_regression_at_tuned_point(testCase)
td  = testCase.TestData;
scn = sim_scenario(td.scn_cfg_s1, td.aj, td.sc);           % sigma_s_db = 0, the P2 point
gap = steady_state_gap(td, td.aj, scn, td.acfg_adapt);
verifyLessThan(testCase, gap, 1.0, ...
    sprintf('adaptive loading broke the P2 gate at its tuned point: gap %.3f dB', gap));
end


% ── P9 gate: sigma_s_db sweep, no re-tune, no regression, real win ─

function test_adaptive_loading_sigma_s_db_sweep(testCase)
td = testCase.TestData;
sigma_s_db_list  = [-20, -10, 0, 10, 20, 30, 40];
regression_tol_db = 0.3;    % adaptive must never be meaningfully worse than fixed
for i = 1:numel(sigma_s_db_list)
    aj = td.aj; aj.sigma_s_db = sigma_s_db_list(i);
    scn = sim_scenario(td.scn_cfg_s1, aj, td.sc);
    gap_fixed = steady_state_gap(td, aj, scn, td.acfg_fixed);
    gap_adapt = steady_state_gap(td, aj, scn, td.acfg_adapt);
    verifyLessThan(testCase, gap_adapt - gap_fixed, regression_tol_db, ...
        sprintf('sigma_s_db=%d: adaptive (%.3f dB) worse than fixed (%.3f dB) by more than %.1f dB', ...
        sigma_s_db_list(i), gap_adapt, gap_fixed, regression_tol_db));
end
end


function test_adaptive_loading_wins_at_high_sigma_s_db(testCase)
% The actual point of P9: with ONE untuned loading_factor_db, the gap no
% longer needs a per-scenario re-tune as sigma_s_db grows past the jammer's
% eigenvalue — verified as a clear win over the fixed loading, not just a tie.
td = testCase.TestData;

aj30 = td.aj; aj30.sigma_s_db = 30;
scn30 = sim_scenario(td.scn_cfg_s1, aj30, td.sc);
gap_fixed30 = steady_state_gap(td, aj30, scn30, td.acfg_fixed);
gap_adapt30 = steady_state_gap(td, aj30, scn30, td.acfg_adapt);
verifyGreaterThan(testCase, gap_fixed30 - gap_adapt30, 0.5, ...
    sprintf('sigma_s_db=30: adaptive (%.3f dB) not enough better than fixed (%.3f dB)', ...
    gap_adapt30, gap_fixed30));

aj40 = td.aj; aj40.sigma_s_db = 40;
scn40 = sim_scenario(td.scn_cfg_s1, aj40, td.sc);
gap_fixed40 = steady_state_gap(td, aj40, scn40, td.acfg_fixed);
gap_adapt40 = steady_state_gap(td, aj40, scn40, td.acfg_adapt);
verifyGreaterThan(testCase, gap_fixed40 - gap_adapt40, 1.5, ...
    sprintf('sigma_s_db=40: adaptive (%.3f dB) not enough better than fixed (%.3f dB)', ...
    gap_adapt40, gap_fixed40));
end
