function tests = test_metrics
% TEST_METRICS  Verify compute_directivity_dbi_grid, compute_hpbw and
% evaluate_metrics against Python references (evaluate_metrics on real CST data).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
root = fileparts(fileparts(here));
testCase.TestData.dir     = here;
testCase.TestData.dataDir = fullfile(root, 'data', 'spacing0.6');
end


function ref = load_ref(testCase, name)
ref = jsondecode(fileread(fullfile(testCase.TestData.dir, 'fixtures', [name '.json'])));
end


function test_directivity_grid_matches_python(testCase)
ref = load_ref(testCase, 'directivity_grid');
af  = ref.af_real + 1i * ref.af_imag;

dbi_auto = compute_directivity_dbi_grid(af, ref.theta_deg(:), ref.phi_deg(:));
dbi_norm = compute_directivity_dbi_grid(af, ref.theta_deg(:), ref.phi_deg(:), 2.5);

verifyEqual(testCase, dbi_auto, ref.dbi_auto, 'RelTol', 1e-9, 'AbsTol', 1e-12);
verifyEqual(testCase, dbi_norm, ref.dbi_norm, 'RelTol', 1e-9, 'AbsTol', 1e-12);
end


function test_hpbw_matches_python(testCase)
ref = load_ref(testCase, 'hpbw');

w_a = compute_hpbw(ref.a_cut(:), ref.a_peak0);                   % no wrap
w_b = compute_hpbw(ref.b_cut(:), ref.b_peak0, ref.b_opp(:));    % left wrap
w_c = compute_hpbw(ref.c_cut(:), ref.c_peak0, ref.c_opp(:));    % right wrap

verifyEqual(testCase, w_a, ref.a_width, 'RelTol', 1e-12, 'no-wrap');
verifyEqual(testCase, w_b, ref.b_width, 'RelTol', 1e-12, 'left-wrap');
verifyEqual(testCase, w_c, ref.c_width, 'RelTol', 1e-12, 'right-wrap');
end


function test_evaluate_metrics_matches_python(testCase)
ref = load_ref(testCase, 'evaluate_metrics');

patterns = load_element_patterns(testCase.TestData.dataDir);
theta    = patterns(1).theta_deg;
phi      = patterns(1).phi_deg;
n        = numel(patterns);
[nt, np_] = size(patterns(1).E_abs);

ep = zeros(n, nt, np_);
for i = 1:n
    ep(i, :, :) = get_component(patterns(i), 'Copol').complex;
end
w = ones(n, 1);   % uniform weights, mirrors the Python fixture

directives = { ...
    struct('type','peak','theta',0,'phi',0,'theta_width',20,'phi_width',20,'weight',1.0), ...
    struct('type','null','theta',40,'phi',0,'theta_width',10,'phi_width',20,'weight',2.0) };

m = evaluate_metrics(ep, theta, phi, w, directives, []);

verifyEqual(testCase, m.global_peak_dbi,       ref.global_peak_dbi,       'RelTol', 1e-9);
verifyEqual(testCase, m.global_peak_theta_deg, ref.global_peak_theta_deg, 'AbsTol', 1e-9);
verifyEqual(testCase, m.global_peak_phi_deg,   ref.global_peak_phi_deg,   'AbsTol', 1e-9);
verifyEqual(testCase, m.hpbw_theta_deg,        ref.hpbw_theta_deg,        'RelTol', 1e-9);
verifyEqual(testCase, m.hpbw_phi_deg,          ref.hpbw_phi_deg,          'RelTol', 1e-9);
verifyEqual(testCase, m.total_cost,            ref.total_cost,            'RelTol', 1e-8, 'AbsTol', 1e-14);
verifyEqual(testCase, m.peak_to_null_ratio_db, ref.peak_to_null_ratio_db, 'RelTol', 1e-9);
verifyEqual(testCase, m.iteration_count, 0);

% Per-directive metrics.
verifyEqual(testCase, numel(m.directive_metrics), numel(ref.directive_metrics));
% Peak directive.
verifyEqual(testCase, m.directive_metrics(1).type, 'peak');
verifyEqual(testCase, m.directive_metrics(1).gain_dbi, ref.directive_metrics(1).gain_dbi, 'RelTol', 1e-9);
verifyTrue(testCase, isempty(m.directive_metrics(1).null_depth_db));   % None for peak
verifyEqual(testCase, m.directive_metrics(1).cost_term, ref.directive_metrics(1).cost_term, 'RelTol', 1e-8, 'AbsTol', 1e-14);
% Null directive.
verifyEqual(testCase, m.directive_metrics(2).type, 'null');
verifyEqual(testCase, m.directive_metrics(2).gain_dbi, ref.directive_metrics(2).gain_dbi, 'RelTol', 1e-9);
verifyEqual(testCase, m.directive_metrics(2).null_depth_db, ref.directive_metrics(2).null_depth_db, 'RelTol', 1e-9);
verifyEqual(testCase, m.directive_metrics(2).cost_term, ref.directive_metrics(2).cost_term, 'RelTol', 1e-8, 'AbsTol', 1e-14);
end
