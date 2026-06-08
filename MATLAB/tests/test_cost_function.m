function tests = test_cost_function
% TEST_COST_FUNCTION  Verify build_directive_physical_masks and build_cost_function
% against Python references (incl. pole-crossing, phi wrap, all aggregations,
% and standard/phase_only/amplitude_only/total modes).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
testCase.TestData.dir = here;
end


function directives = sample_directives()
% Mirror of _DIRECTIVES in gen_reference_fixtures.py (keep in sync).
directives = { ...
    struct('type','peak','theta',30,'phi',90,'theta_width',20,'phi_width',60,'weight',1.0), ...
    struct('type','null','theta',-30,'phi',0,'theta_width',10,'phi_width',20,'weight',3.0,'aggregation','max'), ...
    struct('type','peak','theta',0,'phi',0,'theta_width',10,'phi_width',10,'weight',2.0,'aggregation','min'), ...
    struct('type','null','theta',90,'phi',350,'theta_width',20,'phi_width',40,'weight',5.0), ...
    struct('type','null','theta',185,'phi',200,'theta_width',20,'phi_width',40,'weight',1.0) };
end


function ref = load_ref(testCase, name)
ref = jsondecode(fileread(fullfile(testCase.TestData.dir, 'fixtures', [name '.json'])));
end


function test_physical_masks_match_python(testCase)
ref   = load_ref(testCase, 'phys_masks');
theta = ref.theta_deg(:);
phi   = ref.phi_deg(:);

masks = build_directive_physical_masks(theta, phi, sample_directives());
verifyEqual(testCase, numel(masks), numel(ref.true_counts));

for k = 1:numel(masks)
    expected = squeeze(ref.masks(k, :, :)) == 1;
    verifySize(testCase, masks{k}, [numel(theta), numel(phi)]);
    verifyEqual(testCase, masks{k}, expected, sprintf('mask %d', k));
    verifyEqual(testCase, nnz(masks{k}), ref.true_counts(k), sprintf('count %d', k));
end
end


function test_cost_values_match_python(testCase)
ref   = load_ref(testCase, 'cost_function');
theta = ref.theta_deg(:);
phi   = ref.phi_deg(:);
ep    = ref.ep_real  + 1i * ref.ep_imag;    % (N, N_theta, N_phi)
ep2   = ref.ep2_real + 1i * ref.ep2_imag;
d     = sample_directives();

cf_standard = build_cost_function(ep, theta, phi, d, 'standard');
cf_phase    = build_cost_function(ep, theta, phi, d, 'phase_only');
cf_amp      = build_cost_function(ep, theta, phi, d, 'amplitude_only');
cf_total    = build_cost_function(ep, theta, phi, d, 'standard', ep2);

j_standard = cf_standard(ref.x(:));
j_phase    = cf_phase(ref.x(:));
j_amp      = cf_amp(ref.x_amp(:));
j_total    = cf_total(ref.x(:));

verifyEqual(testCase, j_standard, ref.j_standard, 'RelTol', 1e-9, 'standard');
verifyEqual(testCase, j_phase,    ref.j_phase,    'RelTol', 1e-9, 'phase_only');
verifyEqual(testCase, j_amp,      ref.j_amp,      'RelTol', 1e-9, 'amplitude_only');
verifyEqual(testCase, j_total,    ref.j_total,    'RelTol', 1e-9, 'total');
end


function test_invalid_mode_errors(testCase)
d = sample_directives();
ep = zeros(2, 19, 12);
verifyError(testCase, @() build_cost_function(ep, (0:10:180)', (0:30:330)', d, 'bogus'), ...
    'build_cost_function:BadMode');
end
