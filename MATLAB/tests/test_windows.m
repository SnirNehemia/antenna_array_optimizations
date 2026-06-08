function tests = test_windows
% TEST_WINDOWS  Verify window_1d matches numpy/scipy windows for several lengths.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
testCase.TestData.ref = jsondecode(fileread(fullfile(here, 'fixtures', 'windows.json')));
end


function test_all_windows_match_python(testCase)
ref   = testCase.TestData.ref;
names = cellstr(ref.names);
for ci = 1:numel(ref.cases)
    n = ref.cases(ci).n;
    for wi = 1:numel(names)
        name = names{wi};
        w = window_1d(name, n);
        expected = ref.cases(ci).windows.(name);
        verifySize(testCase, w, [n, 1]);
        verifyEqual(testCase, w, expected(:), 'RelTol', 1e-9, 'AbsTol', 1e-12, ...
            sprintf('%s n=%d', name, n));
    end
end
end


function test_uniform_and_orientation(testCase)
% uniform is all ones; output is always a column vector.
w = window_1d('uniform', 5);
verifyEqual(testCase, w, ones(5, 1));
verifyError(testCase, @() window_1d('bogus', 4), 'window_1d:UnknownWindow');
end
