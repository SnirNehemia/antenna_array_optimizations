function tests = test_weights_x
% TEST_WEIGHTS_X  Unit tests for x_to_weights / weights_to_x against Python refs.
%
%   Run via:  runtests('test_weights_x')  (or the MATLAB MCP test tool).
%   Compares the MATLAB port to fixtures/weights_x.json produced by the Python
%   src/cost/cost_function.py implementation.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
% Put the MATLAB/ port folder (parent of tests/) on the path and load the fixture.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
fixture = fullfile(here, 'fixtures', 'weights_x.json');
testCase.TestData.ref = jsondecode(fileread(fixture));
end


function test_weights_to_x_matches_python(testCase)
ref = testCase.TestData.ref;
w   = ref.w_real(:) + 1i * ref.w_imag(:);

x = weights_to_x(w);

verifySize(testCase, x, [2 * numel(w), 1]);
verifyEqual(testCase, x, ref.x(:), 'AbsTol', 1e-12);
end


function test_x_to_weights_matches_python(testCase)
ref = testCase.TestData.ref;

w = x_to_weights(ref.x_in(:));

verifySize(testCase, w, [numel(ref.x_in) / 2, 1]);
verifyEqual(testCase, real(w), ref.w_decode_real(:), 'AbsTol', 1e-12);
verifyEqual(testCase, imag(w), ref.w_decode_imag(:), 'AbsTol', 1e-12);
end


function test_round_trip_identity(testCase)
% x_to_weights(weights_to_x(w)) == w for an arbitrary complex vector.
w = [3 - 1i; 0 + 0i; -2.5 + 7.25i; 1.0 + 0.5i];
verifyEqual(testCase, x_to_weights(weights_to_x(w)), w, 'AbsTol', 1e-12);
end
