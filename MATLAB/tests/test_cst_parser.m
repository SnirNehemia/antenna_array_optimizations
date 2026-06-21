function tests = test_cst_parser
% TEST_CST_PARSER  Validate parse_cst_file / load_element_patterns against
% Python references on real CST data (data/spacing0.9).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
root = fileparts(fileparts(here));                 % repo root
testCase.TestData.dataDir  = fullfile(root, 'data', 'spacing0.9');
testCase.TestData.oneFile  = fullfile(testCase.TestData.dataDir, ...
                                      'farfield (f=2400) [1].txt');
testCase.TestData.ref = jsondecode(fileread(fullfile(here, 'fixtures', 'cst_parser.json')));
end


function test_grid_shape_and_axes(testCase)
ref = testCase.TestData.ref;
p   = parse_cst_file(testCase.TestData.oneFile);

% Validate parser output (Hard Rule #3): shapes + grid vectors.
fprintf('\n[parser] E_abs size = %s ; n_theta=%d n_phi=%d\n', ...
        mat2str(size(p.E_abs)), ref.n_theta, ref.n_phi);

verifySize(testCase, p.E_abs, [ref.n_theta, ref.n_phi]);
verifyEqual(testCase, numel(p.theta_deg), ref.n_theta);
verifyEqual(testCase, numel(p.phi_deg),   ref.n_phi);
verifyEqual(testCase, p.theta_deg, ref.theta_deg(:), 'AbsTol', 1e-9);
verifyEqual(testCase, p.phi_deg,   ref.phi_deg(:),   'AbsTol', 1e-9);
verifyEqual(testCase, sort(fieldnames(p.components)), sort(ref.component_names(:)));
end


function test_sample_points_match_python(testCase)
ref = testCase.TestData.ref;
p   = parse_cst_file(testCase.TestData.oneFile);
copol = get_component(p, 'Copol');
cross = get_component(p, 'Cross');

for i = 1:numel(ref.samples)
    s = ref.samples(i);
    t = s.t + 1;   % 0-based -> 1-based
    q = s.q + 1;
    msg = sprintf('sample (t=%d,q=%d)', s.t, s.q);
    verifyEqual(testCase, p.E_abs(t, q),       s.E_abs,       'RelTol', 1e-9, msg);
    verifyEqual(testCase, copol.abs(t, q),     s.copol_abs,   'RelTol', 1e-9, msg);
    verifyEqual(testCase, copol.phase(t, q),   s.copol_phase, 'RelTol', 1e-9, msg);
    verifyEqual(testCase, cross.abs(t, q),     s.cross_abs,   'RelTol', 1e-9, msg);
    verifyEqual(testCase, cross.phase(t, q),   s.cross_phase, 'RelTol', 1e-9, msg);
    verifyEqual(testCase, real(copol.complex(t, q)), s.Ecx_real, 'AbsTol', 1e-9, msg);
    verifyEqual(testCase, imag(copol.complex(t, q)), s.Ecx_imag, 'AbsTol', 1e-9, msg);
    verifyEqual(testCase, real(cross.complex(t, q)), s.Xcx_real, 'AbsTol', 1e-9, msg);
    verifyEqual(testCase, imag(cross.complex(t, q)), s.Xcx_imag, 'AbsTol', 1e-9, msg);
end
end


function test_full_slices_match_python(testCase)
% End-to-end check of both axis orderings via a full column and a full row.
ref = testCase.TestData.ref;
p   = parse_cst_file(testCase.TestData.oneFile);
copol = get_component(p, 'Copol');

col0 = copol.complex(:, 1);     % all theta at phi index 1
row0 = copol.complex(1, :).';   % all phi at theta index 1 (as column)

verifyEqual(testCase, real(col0), ref.col0_real(:), 'AbsTol', 1e-9);
verifyEqual(testCase, imag(col0), ref.col0_imag(:), 'AbsTol', 1e-9);
verifyEqual(testCase, real(row0), ref.row0_real(:), 'AbsTol', 1e-9);
verifyEqual(testCase, imag(row0), ref.row0_imag(:), 'AbsTol', 1e-9);
end


function test_load_element_patterns(testCase)
ref      = testCase.TestData.ref;
patterns = load_element_patterns(testCase.TestData.dataDir);

verifyEqual(testCase, numel(patterns), ref.n_elements);

% Each element shares the same grid shape.
sz = size(patterns(1).E_abs);
for i = 2:numel(patterns)
    verifyEqual(testCase, size(patterns(i).E_abs), sz);
end

% Element [1] (first by index) matches the standalone parse of that file.
p1 = parse_cst_file(testCase.TestData.oneFile);
copol1 = get_component(patterns(1), 'Copol');
copol_p1 = get_component(p1, 'Copol');
verifyEqual(testCase, sum(copol1.abs(:)), ref.copol_abs_sum, 'RelTol', 1e-9);
verifyEqual(testCase, copol1.complex, copol_p1.complex, 'AbsTol', 1e-12);
end
