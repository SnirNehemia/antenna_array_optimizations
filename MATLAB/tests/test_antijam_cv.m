function tests = test_antijam_cv
% TEST_ANTIJAM_CV  [P12b] gates for the constant-velocity DoA predictor.
%
%   1. Contract: adapt_cv_init rejects a cv block with any key missing
%      (CLAUDE.md rule 4 — no silent defaults).
%   2. Noise-free tracking: on a clean constant-velocity angle sequence the
%      filter recovers the drift rate, and its lead_steps-ahead prediction
%      lands on the true FUTURE angle rather than the current one.
%   3. The static guard: with a stationary source the estimated speed stays
%      below min_speed_deg_s, so pred.valid is false and a caller keeps its
%      existing branch. This is what makes the CV path a strict addition.
%   4. Outlier gate + track drop: a measurement beyond gate_deg is rejected
%      rather than fused, and a run of misses past max_misses drops the track.
%   5. Mirror fold: on a theta-mirror-degenerate array the filter detects the
%      degeneracy and folds measurements to the branch nearer its prediction,
%      so a branch-hopping MUSIC sequence still yields the correct velocity.
%      On a non-degenerate array the fold stays off.
%   6. Azimuth wrap: a jammer drifting through phi = 0 is tracked without the
%      360 deg residual that an unwrapped innovation would produce.
%
% The filter is exercised directly (not through closed_loop_run) so a failure
% localizes to the estimator rather than to the loop around it.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));

td = struct();
td.dt_s = 0.05;
td.cv = struct('lead_steps', 14, 'q_accel_deg_s2', 0.5, 'r_meas_deg', 1.0, ...
    'gate_deg', 10.0, 'min_track_steps', 20, 'min_speed_deg_s', 0.5, ...
    'max_misses', 20);

% Two toy far-field stacks on a coarse grid. Only the MIRROR SYMMETRY of the
% stack matters to adapt_cv_init — it reads nothing else from them — so these
% are built to have exactly the property under test and nothing more.
td.theta_deg = 0:5:180;
td.phi_deg   = 0:5:355;
n_el = 8;
n_th = numel(td.theta_deg);
n_ph = numel(td.phi_deg);

% (a) Mirror-degenerate: values depend on theta only through |theta - 90|, so
%     e(theta, phi) == e(180 - theta, phi) exactly, like a planar array in the
%     z = 0 plane with no geometric phase re-added.
Em = zeros(n_el, n_th, n_ph);
for it = 1:n_th
    for ip = 1:n_ph
        u = abs(td.theta_deg(it) - 90) / 90;
        v = td.phi_deg(ip) / 360;
        Em(:, it, ip) = exp(1i * 2 * pi * (0:n_el - 1)' * (0.5 * u + 0.3 * v));
    end
end
td.E_mirror = Em;

% (b) Non-degenerate: a genuine theta dependence breaks the symmetry.
Ea = zeros(n_el, n_th, n_ph);
for it = 1:n_th
    for ip = 1:n_ph
        u = td.theta_deg(it) / 180;
        v = td.phi_deg(ip) / 360;
        Ea(:, it, ip) = exp(1i * 2 * pi * (0:n_el - 1)' * (0.7 * u + 0.3 * v));
    end
end
td.E_plain = Ea;

testCase.TestData = td;
end


function test_missing_key_raises(testCase)
% Gate 1 — every cv key is required; none may be silently defaulted.
td = testCase.TestData;
keys = fieldnames(td.cv);
for i = 1:numel(keys)
    bad = rmfield(td.cv, keys{i});
    f = @() adapt_cv_init(bad, td.dt_s, td.E_plain, [], td.theta_deg, td.phi_deg);
    verifyError(testCase, f, 'adapt_cv_init:MissingKey', ...
        sprintf('Dropping cv key ''%s'' should raise.', keys{i}));
end
verifyError(testCase, ...
    @() adapt_cv_init(td.cv, 0, td.E_plain, [], td.theta_deg, td.phi_deg), ...
    'adapt_cv_init:BadStep');
end


function test_tracks_constant_velocity_and_leads(testCase)
% Gate 2 — recover the rate, and predict the FUTURE angle, not the current one.
td = testCase.TestData;
st = adapt_cv_init(td.cv, td.dt_s, td.E_plain, [], td.theta_deg, td.phi_deg);

rate_deg_s = 2.0;                                  % the milestone's drift rate
th0 = 100.0;
n = 200;
for k = 1:n
    t = (k - 1) * td.dt_s;
    [pred, st] = adapt_cv_update(st, th0 + rate_deg_s * t, 260.0, true);
end

verifyTrue(testCase, pred.valid, 'Track should be valid after 200 clean updates.');
verifyEqual(testCase, pred.speed_deg_s, rate_deg_s, 'RelTol', 0.05, ...
    'Estimated angular speed should match the true drift rate within 5%.');

% The whole point of the lead: the prediction must sit lead_steps ahead of the
% current truth, not on it.
t_now    = (n - 1) * td.dt_s;
th_now   = th0 + rate_deg_s * t_now;
th_ahead = th_now + rate_deg_s * td.cv.lead_steps * td.dt_s;
verifyEqual(testCase, pred.theta_deg, th_ahead, 'AbsTol', 0.15, ...
    'Prediction should land on the angle lead_steps into the future.');
verifyEqual(testCase, pred.theta_now_deg, th_now, 'AbsTol', 0.15, ...
    'The filtered CURRENT angle should still match truth now.');
end


function test_static_guard_keeps_cv_out(testCase)
% Gate 3 — a stationary source must never set pred.valid, which is what makes
% STATIC/WINDOW runs bit-identical with the CV block on.
td = testCase.TestData;
st = adapt_cv_init(td.cv, td.dt_s, td.E_plain, [], td.theta_deg, td.phi_deg);
rng(4242, 'twister');
any_valid = false;
for k = 1:300
    % Stationary, with realistic grid-quantized jitter on the measurement.
    th_m = 120.0 + 1.0 * round(randn * 0.4);
    [pred, st] = adapt_cv_update(st, th_m, 260.0, true);
    any_valid = any_valid || pred.valid;
end
verifyFalse(testCase, any_valid, ...
    ['A stationary jammer must never produce a valid CV steer: the speed ' ...
     'gate is what guarantees the non-drifting scenarios are unchanged.']);
end


function test_outlier_gate_and_track_drop(testCase)
% Gate 4 — a wild measurement is rejected, and a long miss run drops the track.
td = testCase.TestData;
st = adapt_cv_init(td.cv, td.dt_s, td.E_plain, [], td.theta_deg, td.phi_deg);
rate = 2.0; th0 = 100.0;
for k = 1:100
    [~, st] = adapt_cv_update(st, th0 + rate * (k - 1) * td.dt_s, 260.0, true);
end
th_before = st.x(1);

% A measurement far outside the gate must not move the state.
[pred, st] = adapt_cv_update(st, th0 + 90.0, 260.0, true);
verifyFalse(testCase, pred.accepted, 'A measurement beyond gate_deg must be rejected.');
verifyEqual(testCase, st.x(1), th_before + st.x(2) * td.dt_s, 'AbsTol', 1e-9, ...
    'A gated-out measurement must leave the state on its predicted trajectory.');

% max_misses consecutive absences drop the track entirely.
for k = 1:(td.cv.max_misses + 2)
    [pred, st] = adapt_cv_update(st, NaN, NaN, false);
end
verifyFalse(testCase, pred.valid, 'Track must be dropped after max_misses.');
verifyEmpty(testCase, st.x, 'A dropped track must clear its state.');
end


function test_mirror_fold(testCase)
% Gate 5 — the degeneracy is detected, and a branch-hopping measurement
% sequence still yields the correct velocity.
td = testCase.TestData;

st_plain = adapt_cv_init(td.cv, td.dt_s, td.E_plain, [], td.theta_deg, td.phi_deg);
verifyFalse(testCase, st_plain.mirror_fold, ...
    'A non-degenerate array must not enable the mirror fold.');

st = adapt_cv_init(td.cv, td.dt_s, td.E_mirror, [], td.theta_deg, td.phi_deg);
verifyTrue(testCase, st.mirror_fold, ...
    'An exactly theta-symmetric array must enable the mirror fold.');
verifyGreaterThan(testCase, st.mirror_coh, 0.99);

% Truth drifts 110 -> 130 deg; MUSIC reports the MIRROR branch on alternate
% steps, exactly as it does on the real degenerate array.
rate = 2.0; th0 = 110.0;
for k = 1:200
    th_true = th0 + rate * (k - 1) * td.dt_s;
    th_meas = th_true;
    if mod(k, 2) == 0
        th_meas = 180 - th_true;               % the ambiguous branch
    end
    [pred, st] = adapt_cv_update(st, th_meas, 260.0, true);
end
verifyTrue(testCase, pred.valid);
verifyEqual(testCase, pred.speed_deg_s, rate, 'RelTol', 0.10, ...
    ['With the fold on, a branch-hopping DoA sequence must still give the ' ...
     'correct drift rate — that is the whole reason the fold exists.']);
end


function test_azimuth_wrap(testCase)
% Gate 6 — drifting through phi = 0 must not trip the gate on a 359 deg
% residual, and the reported azimuth must stay inside [0, 360).
td = testCase.TestData;
st = adapt_cv_init(td.cv, td.dt_s, td.E_plain, [], td.theta_deg, td.phi_deg);
rate = 2.0;
ph0  = 355.0;
for k = 1:200
    ph_true = mod(ph0 + rate * (k - 1) * td.dt_s, 360);
    [pred, st] = adapt_cv_update(st, 90.0, ph_true, true);
end
verifyTrue(testCase, pred.valid, 'Track must survive the phi = 0 crossing.');
verifyEqual(testCase, pred.speed_deg_s, rate, 'RelTol', 0.10, ...
    'Azimuth drift rate must be recovered across the wrap.');
verifyGreaterThanOrEqual(testCase, pred.phi_deg, 0);
verifyLessThan(testCase, pred.phi_deg, 360);
end
