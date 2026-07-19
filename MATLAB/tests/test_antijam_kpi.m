function tests = test_antijam_kpi
% TEST_ANTIJAM_KPI  Unit tests for kpi_evaluate on hand-computable logs
% (P6). Uses a 2-element toy cut where every KPI has a closed-form value.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));
end


function [run_log, scenario, aj] = toy_case()
% 2-el ULA cut over -90..90 (61 points); jammer static at +40 deg, one
% turn_on event at t = 1 s; T = 5 steps of dt = 1 s.
ang = linspace(-90, 90, 61);
E   = exp(1i * pi * (0:1).' * sind(ang));       % 2-el half-wavelength ULA
e_s = E(:, 31);                                  % boresight
T   = 5;
scenario = struct('id', 'toy', 't_s', 0:T-1, ...
    'theta_j_deg', 40 * ones(1, T), ...
    'jammer_on',  logical([0 1 1 1 1]), ...
    'jn_ratio_db', 20 * ones(1, T), ...
    'events', {{struct('t_s', 1.0, 'type', 'turn_on')}});
aj = struct('sinr_min_db', 5.0, 'theta_s_deg', 0.0);

% Weights: steps 1-2 quiescent (matched filter), steps 3-5 = MVDR against
% the jammer -> exact null at +40, SINR above threshold.
w_q = adapt_lcmv(eye(2), e_s, 0);
R   = 100 * (E(:, nearest_index(ang(:), 40)) * E(:, nearest_index(ang(:), 40))') + eye(2);
w_n = adapt_lcmv(R, e_s, 0);
W   = [w_q, w_q, w_n, w_n, w_n];

run_log = struct();
run_log.W = W;
run_log.cut = struct('E1', E, 'E2', [], 'angle_deg', ang, 'e_s', e_s);
% SINR timeline chosen by hand around the 5 dB threshold:
run_log.sinr_db        = [10, -3, 2, 7, 8];     % dips at turn-on, recovers @ k=4
run_log.oracle_sinr_db = [10,  8, 8, 8, 8];
end


function test_availability_recovery_gap(testCase)
[run_log, scenario, aj] = toy_case();
kpi = kpi_evaluate(run_log, scenario, aj);
% 3 of 5 steps >= 5 dB.
verifyEqual(testCase, kpi.availability, 3 / 5, 'AbsTol', 1e-12);
% Event at t=1 (k=2): sinr [-3, 2, 7, ...] -> first >= thr at offset 2.
verifyEqual(testCase, kpi.recovery_steps, 2);
% Oracle gap element-wise + mean.
verifyEqual(testCase, kpi.oracle_gap_db, [0, 11, 6, 1, 0], 'AbsTol', 1e-12);
verifyEqual(testCase, kpi.oracle_gap_mean_db, mean([0, 11, 6, 1, 0]), 'AbsTol', 1e-12);
end


function test_null_pointing_and_gain_penalty(testCase)
[run_log, scenario, aj] = toy_case();
kpi = kpi_evaluate(run_log, scenario, aj);
% Jammer off at k=1 -> NaN. MVDR steps: nearest pattern null within one
% grid cell (3 deg) of +40. (The quiescent 2-el pattern also has SOME null,
% so k=2 is finite but far-off is fine — only check the MVDR steps.)
verifyTrue(testCase, isnan(kpi.null_pointing_err_deg(1)));
verifyLessThanOrEqual(testCase, max(kpi.null_pointing_err_deg(3:5)), 3.0);
% Quiescent steps have zero penalty by construction (w = reference).
verifyEqual(testCase, kpi.peak_gain_penalty_db(1:2), zeros(1, 2), 'AbsTol', 1e-9);
% MVDR-with-jammer steps sacrifice a little gain, never gain any.
verifyLessThanOrEqual(testCase, kpi.peak_gain_penalty_db(3:5), zeros(1, 3) + 1e-9);
% Missing threshold key raises.
verifyError(testCase, @() kpi_evaluate(run_log, scenario, struct('theta_s_deg', 0)), ...
    'kpi_evaluate:MissingKey');
end
