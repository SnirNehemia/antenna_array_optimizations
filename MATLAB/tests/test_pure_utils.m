function tests = test_pure_utils
% TEST_PURE_UTILS  Unit tests for nearest_index, angular_window_mask,
% and get_directive_field against Python references where applicable.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
testCase.TestData.dir = here;
end


function ref = load_ref(testCase, name)
fixture = fullfile(testCase.TestData.dir, 'fixtures', [name '.json']);
ref = jsondecode(fileread(fixture));
end


function test_nearest_index_matches_python(testCase)
ref = load_ref(testCase, 'nearest_index');
grid = ref.grid_deg(:);
for i = 1:numel(ref.targets)
    idx = nearest_index(grid, ref.targets(i));
    % MATLAB is 1-based; Python argmin is 0-based.
    verifyEqual(testCase, idx, ref.idx0_python(i) + 1, ...
        sprintf('target %.2f', ref.targets(i)));
end
end


function test_angular_window_mask_matches_python(testCase)
ref   = load_ref(testCase, 'angular_window_mask');
theta = ref.theta_deg(:);
phi   = ref.phi_deg(:);
for k = 1:numel(ref.cases)
    c = ref.cases(k);
    mask = angular_window_mask(theta, phi, c.tt, c.tp, c.tw, c.pw);
    expected = squeeze(ref.masks(k, :, :)) == 1;   % (N_theta x N_phi) logical
    verifySize(testCase, mask, [numel(theta), numel(phi)]);
    verifyEqual(testCase, mask, expected, sprintf('case %d', k));
end
end


function test_get_directive_field(testCase)
% Mirrors Python dict.get: present field returns value, absent returns default.
d = struct('type', 'peak', 'theta', 30.0, 'weight', 5.0);
verifyEqual(testCase, get_directive_field(d, 'theta', -1), 30.0);
verifyEqual(testCase, get_directive_field(d, 'weight', 1.0), 5.0);
verifyEqual(testCase, get_directive_field(d, 'phi', 0.0), 0.0);        % absent -> default
verifyEqual(testCase, get_directive_field(d, 'theta_width', 7.5), 7.5); % absent -> default
end
