function tests = test_compare_classical
% TEST_COMPARE_CLASSICAL  Verify the compare_classical math helpers and
% run_scenario (classical techniques) against Python references.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
testCase.TestData.dir = here;
end


function ref = load_ref(testCase, name)
ref = jsondecode(fileread(fullfile(testCase.TestData.dir, 'fixtures', [name '.json'])));
end


function test_helpers_match_python(testCase)
ref   = load_ref(testCase, 'compare_helpers');
theta = ref.theta_deg(:);
phi   = ref.phi_deg(:);
n_side = 3; d = 0.5;
st = ref.steer_theta; sp = ref.steer_phi;

% build_ura_element_patterns
ep = build_ura_element_patterns(n_side, d, theta, phi);
ep_ref = ref.ep_real + 1i * ref.ep_imag;
verifyEqual(testCase, ep, ep_ref, 'AbsTol', 1e-12, 'build_ura');

% steering vectors
sgeo = steering_phase_vector(n_side, d, st, sp);
verifyEqual(testCase, sgeo, ref.steering_geo_real(:) + 1i * ref.steering_geo_imag(:), ...
    'AbsTol', 1e-12, 'steering_geo');
sdd = data_driven_steering_vector(ep, theta, phi, st, sp);
verifyEqual(testCase, sdd, ref.steering_dd_real(:) + 1i * ref.steering_dd_imag(:), ...
    'AbsTol', 1e-12, 'steering_dd');

% classical_weights (data-driven)
names = cellstr(ref.cw_names);
for i = 1:numel(names)
    nm = names{i};
    w  = classical_weights(nm, n_side, d, st, sp, ep, theta, phi);
    expected = ref.cw.(nm).real(:) + 1i * ref.cw.(nm).imag(:);
    verifyEqual(testCase, w, expected, 'AbsTol', 1e-12, nm);
end

% principal_plane_theta_axis
verifyEqual(testCase, principal_plane_theta_axis(theta), ref.theta_axis(:), 'AbsTol', 1e-12);

% principal_plane_cut + evaluate_directive_metrics on a concrete pattern.
w_ham  = classical_weights('hamming', n_side, d, st, sp, ep, theta, phi);
w_norm = power_normalize_weights(w_ham);
af     = compute_array_factor(w_norm, ep);
af_mag = abs(af);
verifyEqual(testCase, af_mag, ref.af_mag, 'AbsTol', 1e-12, 'af_mag');
verifyEqual(testCase, principal_plane_cut(af, phi), ref.cut(:), 'AbsTol', 1e-12, 'cut');

edm = evaluate_directive_metrics(af_mag, theta, phi, ref.directive);
verifyEqual(testCase, edm.peak_dbi, ref.edm.peak_dbi, 'RelTol', 1e-9, 'AbsTol', 1e-12);
verifyEqual(testCase, edm.min_db,   ref.edm.min_db,   'RelTol', 1e-9, 'AbsTol', 1e-12);
verifyEqual(testCase, edm.mean_db,  ref.edm.mean_db,  'RelTol', 1e-9, 'AbsTol', 1e-12);
verifyEqual(testCase, edm.max_db,   ref.edm.max_db,   'RelTol', 1e-9, 'AbsTol', 1e-12);

% power_normalize_weights
wtest = ref.wtest_real(:) + 1i * ref.wtest_imag(:);
verifyEqual(testCase, power_normalize_weights(wtest), ...
    ref.wtest_norm_real(:) + 1i * ref.wtest_norm_imag(:), 'AbsTol', 1e-12);
end


function test_run_scenario_classical_match_python(testCase)
ref   = load_ref(testCase, 'run_scenario');
theta = (0:10:180)';
phi   = (0:30:330)';
n_side = 3; d = 0.5;
ep = build_ura_element_patterns(n_side, d, theta, phi);

scenario = struct( ...
    'steer_theta_deg', 0.0, 'steer_phi_deg', 0.0, ...
    'directives', {{ ...
        struct('type','peak','theta',0,'phi',0,'theta_width',20,'phi_width',360,'weight',1.0), ...
        struct('type','null','theta',40,'phi',0,'theta_width',10,'phi_width',360,'weight',2.0) }});
cfg = struct('n_restarts', 2, 'max_iterations', 100, 'cost_tolerance', 1e-9);

results = run_scenario(scenario, ep, theta, phi, n_side, d, cfg);

verifyTrue(testCase, isfield(results, 'optimized'));
names = cellstr(ref.names);
for i = 1:numel(names)
    nm = names{i};
    r  = results.(nm);
    verifyEqual(testCase, numel(r.cut_vm), ref.n_cut, sprintf('%s cut len', nm));
    verifyEqual(testCase, numel(r.directive_metrics), ref.n_directives);
    for di = 1:ref.n_directives
        exp_dm = ref.classical.(nm).dm(di);
        got    = r.directive_metrics(di);
        verifyEqual(testCase, got.peak_dbi, exp_dm.peak_dbi, 'RelTol', 1e-9, 'AbsTol', 1e-12);
        verifyEqual(testCase, got.mean_db,  exp_dm.mean_db,  'RelTol', 1e-9, 'AbsTol', 1e-12);
        verifyEqual(testCase, got.min_db,   exp_dm.min_db,   'RelTol', 1e-9, 'AbsTol', 1e-12);
        verifyEqual(testCase, got.max_db,   exp_dm.max_db,   'RelTol', 1e-9, 'AbsTol', 1e-12);
    end
end
% Optimized cut has the same length.
verifyEqual(testCase, numel(results.optimized.cut_vm), ref.n_cut);
end
