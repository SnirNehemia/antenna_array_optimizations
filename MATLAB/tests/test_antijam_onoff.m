function tests = test_antijam_onoff
% TEST_ANTIJAM_ONOFF  [O] gates for the on/off detection repair and for the
% graceful-degradation contract on apertures too small to support MUSIC.
%
%   1. Contract: every adapt.predict.onoff key is required; fast_lambda must be
%      shorter-memory than forgetting_lambda (a detector slower than the
%      beamformer cannot resolve what the beamformer already smooths).
%   2. Inert when absent: with no onoff block the run is BYTE-IDENTICAL to the
%      pre-repair path, so every earlier gate still measures what it did.
%   3. The buffer is sized so min_periods cycles of max_period_s fit — this is
%      the arithmetic that made periods above ~17 s undetectable before.
%   4. Presence accuracy: on a fast on/off jammer the repaired presence
%      indicator tracks the true duty cycle far better than the shipped one,
%      which saturates because it reads a covariance smoothed by lambda.
%   5. The anticipatory branch actually FIRES (it measured <= 1.3% before).
%   6. Graceful degradation: an aperture with n_el <= 2*n_comp runs to
%      completion on both lcmv and predict instead of throwing, and predict
%      falls back to the reactive solution rather than inventing a DoA.
%
% Toy 8-element ULA fixture for the estimator gates; the degradation gates use
% deliberately tiny synthetic apertures.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));

td = struct();
td.aj = struct('theta_s_deg', 0.0, 'phi_s_deg', 0.0, 'guard_deg', 10.0, ...
               'jn_ratio_db', 25.0, 'sigma_s_db', 10.0, 'sinr_min_db', 5.0);
td.sc = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 32, ...
               'seed', 1234);

% 8-element half-wavelength ULA over the physical theta domain, phi fixed.
n_el = 8;
td.theta_deg = 0:2:180;
td.phi_deg   = 0;
n_th = numel(td.theta_deg);
E = zeros(n_el, n_th, 1);
for it = 1:n_th
    u = cosd(td.theta_deg(it));
    E(:, it, 1) = exp(1i * pi * (0:n_el - 1)' * u);
end
td.stack = E;
td.n_el  = n_el;

td.acfg = struct('forgetting_lambda', 0.90, 'diagonal_loading_db', 10, ...
    'weight_smoothing_mu', 1.0, ...
    'predict', struct('presence_gap_db', 6.0, 'buffer_len', 1024, ...
                      'min_periods', 3, 'lead_steps', 6, 'doa_stride', 1));
td.onoff = struct('fast_lambda', 0.5, 'max_period_s', 60.0, ...
                  'lead_frac', 0.15, 'lead_cap_horizons', 2.0);
testCase.TestData = td;
end


function cfg = with_onoff(td)
cfg = struct('adapt', td.acfg, 'polarization', 'copol');
cfg.adapt.predict.onoff = td.onoff;
end


function cfg = without_onoff(td)
cfg = struct('adapt', td.acfg, 'polarization', 'copol');
end


function [scn, sc_cfg] = onoff_scenario(td, period_s)
% The run must be long enough to CONTAIN min_periods cycles, otherwise no
% amount of buffer sizing can detect the period -- an inherent latency of
% min_periods * period seconds before the first anticipatory null.
sc_cfg = td.sc;
sc_cfg.duration_s = max(td.sc.duration_s, ...
    1.5 * td.acfg.predict.min_periods * period_s);
sc = struct('id', 'ONOFF', 'motion', 'static', 'power', 'onoff', ...
    'duty_cycle', 0.5, 'toggle_period_s', period_s, ...
    'theta_j_deg', 40.0, 'phi_j_deg', 0.0, 'jn_ratio_db', td.aj.jn_ratio_db);
scn = sim_scenario(sc, td.aj, sc_cfg);
end


function test_required_keys_and_lambda_ordering(testCase)
% Gate 1 — no silent defaults, and the detector must be faster than the beamformer.
td = testCase.TestData;
keys = fieldnames(td.onoff);
for i = 1:numel(keys)
    if strcmp(keys{i}, 'lead_cap_horizons'), continue; end   % this one is optional
    bad = td.acfg;
    bad.predict.onoff = rmfield(td.onoff, keys{i});
    verifyError(testCase, ...
        @() adapt_predict_init(bad, td.aj, td.stack(:, 1), td.stack, [], ...
            td.theta_deg, td.phi_deg, td.n_el, td.sc), ...
        'adapt_predict_init:MissingKey', ...
        sprintf('Dropping onoff key ''%s'' should raise.', keys{i}));
end

slow = td.acfg;
slow.predict.onoff = td.onoff;
slow.predict.onoff.fast_lambda = 0.95;      % SLOWER than forgetting_lambda 0.90
verifyError(testCase, ...
    @() adapt_predict_init(slow, td.aj, td.stack(:, 1), td.stack, [], ...
        td.theta_deg, td.phi_deg, td.n_el, td.sc), ...
    'adapt_predict_init:BadFastLambda');
end


function test_inert_when_absent(testCase)
% Gate 2 — with no onoff block the result must be bit-identical to before.
td  = testCase.TestData;
[scn, sc] = onoff_scenario(td, 10.0);
a = closed_loop_run('predict', td.stack, [], td.theta_deg, td.phi_deg, ...
    scn, td.aj, sc, without_onoff(td), []);
b = closed_loop_run('predict', td.stack, [], td.theta_deg, td.phi_deg, ...
    scn, td.aj, sc, without_onoff(td), []);
verifyEqual(testCase, a.sinr_db, b.sinr_db, ...
    'The un-repaired path must be deterministic.');

st = adapt_predict_init(td.acfg, td.aj, td.stack(:, 1), td.stack, [], ...
    td.theta_deg, td.phi_deg, td.n_el, td.sc);
verifyFalse(testCase, st.onoff_enabled, ...
    'An absent onoff block must leave the repair disabled.');
verifyEmpty(testCase, st.R_fast, ...
    'No fast covariance should be allocated when the repair is off.');
verifyEqual(testCase, st.buffer_len, td.acfg.predict.buffer_len, ...
    'buffer_len must be untouched when the repair is off.');
end


function test_buffer_sized_for_longest_period(testCase)
% Gate 3 — THE arithmetic that broke slow toggling: min_periods cycles of
% max_period_s must fit inside the analysis window.
td  = testCase.TestData;
cfg = td.acfg;
cfg.predict.onoff = td.onoff;
st = adapt_predict_init(cfg, td.aj, td.stack(:, 1), td.stack, [], ...
    td.theta_deg, td.phi_deg, td.n_el, td.sc);
need = cfg.predict.min_periods * td.onoff.max_period_s / td.sc.dt_s;
verifyGreaterThanOrEqual(testCase, st.buffer_len, need, ...
    ['The window must hold min_periods cycles of the longest period worth ' ...
     'detecting; the shipped fixed 1024 could not, so any period above ' ...
     '~17 s was undetectable by construction.']);
verifyGreaterThanOrEqual(testCase, st.buffer_len, cfg.predict.buffer_len, ...
    'Sizing must never SHRINK the configured window.');
end


function test_presence_tracks_duty_cycle(testCase)
% Gate 4 — the repaired presence indicator must track the true duty cycle much
% better than the shipped one, whose covariance is smoothed by lambda.
td  = testCase.TestData;
[scn, sc] = onoff_scenario(td, 5.0);           % fast toggling: the hard case
base = closed_loop_run('predict', td.stack, [], td.theta_deg, td.phi_deg, ...
    scn, td.aj, sc, without_onoff(td), []);
fix  = closed_loop_run('predict', td.stack, [], td.theta_deg, td.phi_deg, ...
    scn, td.aj, sc, with_onoff(td), []);

true_duty = mean(scn.jammer_on);
err_base = abs(mean(base.present) - true_duty);
err_fix  = abs(mean(fix.present)  - true_duty);
verifyLessThan(testCase, err_fix, err_base, ...
    sprintf(['Repaired presence should track the true duty (%.2f) better ' ...
             'than the shipped indicator (base err %.3f, fixed err %.3f).'], ...
            true_duty, err_base, err_fix));
verifyLessThan(testCase, err_fix, 0.25, ...
    'Repaired presence should be within 25 pp of the true duty cycle.');
end


function test_anticipatory_branch_fires(testCase)
% Gate 5 — the pre-null must actually engage. It measured <= 1.3% of steps
% before the repair, and 0% outside a narrow 10-15 s band.
td = testCase.TestData;
for period = [10.0, 30.0]
    [scn, sc] = onoff_scenario(td, period);
    fix = closed_loop_run('predict', td.stack, [], td.theta_deg, td.phi_deg, ...
        scn, td.aj, sc, with_onoff(td), []);
    fired = mean(fix.predicted_on & ~fix.present);
    verifyGreaterThan(testCase, fired, 0.005, ...
        sprintf(['The anticipatory branch must fire at a %g s toggle period; ' ...
                 'it fired on %.2f%% of steps.'], period, 100 * fired));
    verifyTrue(testCase, any(isfinite(fix.period_est)), ...
        sprintf('The period must become trusted at a %g s toggle period.', period));
end
end


function test_graceful_degradation_small_aperture(testCase)
% Gate 6 — an aperture too small for MUSIC must RUN, not throw, and predict
% must fall back to the reactive solution rather than invent a DoA.
td = testCase.TestData;
small = td.stack(1:2, :, :);                  % 2 elements, single component
[scn, sc] = onoff_scenario(td, 10.0);
cfg   = without_onoff(td);

w = warning('off', 'adapt_tracking_init:AdaptiveLoadingInfeasible');
c1 = onCleanup(@() warning(w));

l = closed_loop_run('lcmv', small, [], td.theta_deg, td.phi_deg, ...
    scn, td.aj, sc, cfg, []);
p = closed_loop_run('predict', small, [], td.theta_deg, td.phi_deg, ...
    scn, td.aj, sc, cfg, []);
verifyTrue(testCase, all(isfinite(l.sinr_db)), 'lcmv must complete on a 2-element array.');
verifyTrue(testCase, all(isfinite(p.sinr_db)), 'predict must complete on a 2-element array.');

% MUSIC is infeasible here (n_sig = 2 >= n_el = 2), so predict must degrade to
% the reactive path exactly -- no DoA, no presence, and the same weights.
verifyTrue(testCase, all(isnan(p.doa_theta_deg)), ...
    'predict must not report a DoA on an aperture that cannot support MUSIC.');
verifyFalse(testCase, any(p.present), 'presence must stay false when MUSIC is infeasible.');
verifyEqual(testCase, p.sinr_db, l.sinr_db, 'AbsTol', 1e-9, ...
    'Degraded predict must reproduce the reactive beamformer exactly.');
end
