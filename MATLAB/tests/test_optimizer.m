function tests = test_optimizer
% TEST_OPTIMIZER  Verify run_optimizer (fmincon) on a small synthetic URA:
% structural correctness, improvement over uniform, and convergence to the same
% global optimum as the Python scipy reference (cross-solver tolerance).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
testCase.TestData.ref = jsondecode(fileread(fullfile(here, 'fixtures', 'optimizer.json')));
end


function [ep, theta, phi, directives, config, ref] = setup_problem(testCase)
ref   = testCase.TestData.ref;
ep    = ref.ep_real + 1i * ref.ep_imag;     % (9, 19, 12)
theta = ref.theta_deg(:);
phi   = ref.phi_deg(:);
directives = {ref.directive};
config     = ref.config;
end


function test_structure_and_run_count(testCase)
[ep, theta, phi, directives, config, ref] = setup_problem(testCase);
result = run_optimizer(ep, theta, phi, directives, config);

verifyEqual(testCase, numel(result.all_cost_histories), ref.n_runs);
verifyEqual(testCase, numel(result.all_run_labels),     ref.n_runs);
verifySize(testCase, result.weights_complex, [ref.n_elements, 1]);
verifyTrue(testCase, result.best_run_index >= 1 && result.best_run_index <= ref.n_runs);
verifyTrue(testCase, iscell(result.weights_history));
verifyEqual(testCase, numel(result.weights_history), numel(result.cost_history));
for i = 1:numel(result.weights_history)
    verifySize(testCase, result.weights_history{i}, [ref.n_elements, 1]);
end
end


function test_converged_cost_matches_python(testCase)
[ep, theta, phi, directives, config, ref] = setup_problem(testCase);
result = run_optimizer(ep, theta, phi, directives, config);

% Optimizer must do at least as well as uniform weights.
verifyLessThanOrEqual(testCase, result.result.fun, ref.j_uniform + 1e-9);

% Same global optimum as scipy (different solver -> modest tolerance).
verifyEqual(testCase, result.result.fun, ref.j_star, 'RelTol', 1e-4, 'AbsTol', 1e-4);
end


function test_determinism(testCase)
[ep, theta, phi, directives, config, ~] = setup_problem(testCase);
r1 = run_optimizer(ep, theta, phi, directives, config);
r2 = run_optimizer(ep, theta, phi, directives, config);
% Seeded rng(0) inside run_optimizer makes repeated calls identical.
verifyEqual(testCase, r2.weights_complex, r1.weights_complex, 'AbsTol', 1e-12);
verifyEqual(testCase, r2.result.fun, r1.result.fun, 'AbsTol', 1e-12);
end
