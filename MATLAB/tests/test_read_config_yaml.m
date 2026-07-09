function tests = test_read_config_yaml
% TEST_READ_CONFIG_YAML  Verify read_config_yaml parses config.yaml and
% test_config.yaml into values matching PyYAML (yaml.safe_load).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
root = fileparts(fileparts(here));
testCase.TestData.root = root;
testCase.TestData.ref  = jsondecode(fileread(fullfile(here, 'fixtures', 'config_yaml.json')));
end


function test_config_yaml(testCase)
% Reads a frozen reference config (fixtures/config_yaml_ref.yaml) rather than
% the repo-root config.yaml, which users edit freely. The fixture JSON was
% generated from this same reference, so the two stay in lockstep.
ref = testCase.TestData.ref.config;
here = fileparts(mfilename('fullpath'));
c   = read_config_yaml(fullfile(here, 'fixtures', 'config_yaml_ref.yaml'));

verifyEqual(testCase, c.element_patterns_dir, ref.element_patterns_dir);
verifyEqual(testCase, c.polarization,         ref.polarization);

verifyTrue(testCase, iscell(c.directives));
verifyEqual(testCase, numel(c.directives), ref.n_directives);

d0 = c.directives{1};
verifyEqual(testCase, d0.type,        ref.d0_type);
verifyEqual(testCase, d0.theta,       ref.d0_theta);
verifyEqual(testCase, d0.phi,         ref.d0_phi);
verifyEqual(testCase, d0.theta_width, ref.d0_theta_width);
verifyEqual(testCase, d0.phi_width,   ref.d0_phi_width);
verifyEqual(testCase, d0.weight,      ref.d0_weight);
verifyTrue(testCase,  isempty(d0.level));          % null -> []

d1 = c.directives{2};
verifyEqual(testCase, d1.type,   ref.d1_type);
verifyEqual(testCase, d1.theta,  ref.d1_theta);
verifyEqual(testCase, d1.weight, ref.d1_weight);

verifyEqual(testCase, c.optimizer.max_iterations, ref.opt_max_iter);
verifyEqual(testCase, c.optimizer.cost_tolerance, ref.opt_ctol, 'RelTol', 1e-12);
verifyEqual(testCase, c.optimizer.n_restarts,     ref.opt_n_restarts);
verifyEqual(testCase, c.optimizer.amplitude_bounds(:), ref.opt_bounds(:), 'AbsTol', 1e-12);
verifyEqual(testCase, logical(c.optimizer.phase_only),              logical(ref.opt_phase_only));
verifyEqual(testCase, logical(c.optimizer.amplitude_only),          logical(ref.opt_amp_only));
verifyEqual(testCase, logical(c.optimizer.use_uniform_init),        logical(ref.opt_use_uniform));
verifyEqual(testCase, logical(c.optimizer.use_single_element_init), logical(ref.opt_use_single));

verifyEqual(testCase, c.output.results_dir,           ref.out_results_dir);
verifyEqual(testCase, c.output.plot_cut_type,         ref.out_cut_type);
verifyEqual(testCase, c.output.plot_dynamic_range_db, ref.out_dyn_range);
verifyEqual(testCase, logical(c.output.save_pattern_gif), logical(ref.out_save_gif));
verifyEqual(testCase, c.output.gif_max_frames,        ref.out_gif_frames);
verifyEqual(testCase, logical(c.output.plot_equal_area), logical(ref.out_equal_area));
end


function test_test_config_yaml(testCase)
ref = testCase.TestData.ref.test_config;
c   = read_config_yaml(fullfile(testCase.TestData.root, 'test_config.yaml'));

verifyEqual(testCase, c.n_side,         ref.n_side);
verifyEqual(testCase, c.d_over_lambda,  ref.d_over_lambda, 'RelTol', 1e-12);
verifyEqual(testCase, c.theta_step_deg, ref.theta_step_deg);
verifyEqual(testCase, c.phi_step_deg,   ref.phi_step_deg);
verifyEqual(testCase, c.element_source, ref.element_source);
verifyEqual(testCase, c.element_patterns_dir, ref.element_patterns_dir);
verifyEqual(testCase, c.polarization,   ref.polarization);
verifyEqual(testCase, c.n_restarts,     ref.n_restarts);
verifyEqual(testCase, c.max_iterations, ref.max_iterations);
verifyEqual(testCase, c.cost_tolerance, ref.cost_tolerance, 'RelTol', 1e-12);
verifyEqual(testCase, c.plot_dynamic_range_db, ref.plot_dynamic_range_db);

verifyTrue(testCase, iscell(c.scenarios));
verifyEqual(testCase, numel(c.scenarios), ref.n_scenarios);

s0 = c.scenarios{1};
verifyEqual(testCase, s0.label,           ref.s0_label);
verifyEqual(testCase, s0.steer_theta_deg, ref.s0_steer_theta);
verifyEqual(testCase, s0.null_theta_deg,  ref.s0_null_theta);
verifyTrue(testCase, iscell(s0.directives));
verifyEqual(testCase, numel(s0.directives), ref.s0_n_directives);

verifyEqual(testCase, s0.directives{1}.type,      ref.s0_d0_type);
verifyEqual(testCase, s0.directives{1}.theta,     ref.s0_d0_theta);
verifyEqual(testCase, s0.directives{1}.phi,       ref.s0_d0_phi);
verifyEqual(testCase, s0.directives{1}.phi_width, ref.s0_d0_phi_width);
verifyEqual(testCase, s0.directives{2}.type,      ref.s0_d1_type);
verifyEqual(testCase, s0.directives{2}.theta,     ref.s0_d1_theta);
verifyEqual(testCase, s0.directives{2}.weight,    ref.s0_d1_weight);
end
