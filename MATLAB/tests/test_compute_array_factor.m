function tests = test_compute_array_factor
% TEST_COMPUTE_ARRAY_FACTOR  Verify the array factor against the Python reference.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
fixture = fullfile(here, 'fixtures', 'compute_array_factor.json');
testCase.TestData.ref = jsondecode(fileread(fixture));
end


function test_af_matches_python(testCase)
ref = testCase.TestData.ref;

% Reconstruct the (N_elements, N_theta, N_phi) complex stack and weights.
ep = ref.ep_real + 1i * ref.ep_imag;     % 3D (N, N_theta, N_phi)
w  = ref.w_real(:) + 1i * ref.w_imag(:); % (N x 1)

af = compute_array_factor(w, ep);

expected = ref.af_real + 1i * ref.af_imag;   % (N_theta x N_phi)
verifySize(testCase, af, size(expected));
verifyEqual(testCase, af, expected, 'AbsTol', 1e-12, 'RelTol', 1e-12);
end


function test_uniform_weights_sum(testCase)
% Sanity: with unit weights the AF is the plain element-wise sum over elements.
ref = testCase.TestData.ref;
ep  = ref.ep_real + 1i * ref.ep_imag;
n   = size(ep, 1);
af  = compute_array_factor(ones(n, 1), ep);
verifyEqual(testCase, af, squeeze(sum(ep, 1)), 'AbsTol', 1e-12);
end
