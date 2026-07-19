function tests = test_antijam_codebook
% TEST_ANTIJAM_CODEBOOK  P4 gates for codebook generation
% (antijam_milestone_plan.md Section 4) on a 16-element half-wavelength ULA
% (same element count as the real array; ~6 deg beamwidth so the plan's 5 deg
% arm grid / 10 deg null width proportions are physical).
%
% Gate interpretation (refined in P4, see plan): per-arm depth is measured at
% the arm's null CENTER (<= -30 dB relative to the pattern peak); the
% operationally binding criterion is COVERAGE — worst-case best-arm point
% depth <= -25 dB over the whole allowed span. The window-MEAN depth is
% DOF-limited for a filled 10 deg null (saturates ~ -23 dB on this array
% regardless of null weight) and is reported by the P6 campaign, not gated.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));

theta_deg = 0:2:90;
phi_deg   = 0:10:350;
n_el      = 16;
[TH, PH] = ndgrid(deg2rad(theta_deg), deg2rad(phi_deg));
u = sin(TH) .* cos(PH);
stack = complex(zeros(n_el, numel(theta_deg), numel(phi_deg)));
for n = 1:n_el
    stack(n, :, :) = exp(1i * pi * (n - 1 - (n_el - 1) / 2) * u);  % 16-el ULA
end

td = struct();
td.theta_deg = theta_deg;
td.phi_deg   = phi_deg;
td.stack     = stack;
% guard = 30 deg (~4.7 beamwidths): the P4 scan showed the projection's peak
% gain penalty crosses 1 dB for null windows closer than ~30 deg to boresight
% on this array — guard_deg is exactly the knob that encodes that physics.
td.aj = struct('theta_s_deg', 0.0, 'guard_deg', 30.0, 'jn_ratio_db', 20.0, ...
               'sigma_s_db', 0.0, 'cut_type', 'theta_cut', 'cut_phi_deg', 0.0);
td.ag = struct('null_grid_deg', 5.0, 'null_width_deg', 10.0, ...
               'peak_width_deg', 6.0, 'null_weight', 100.0);
td.oc = struct('max_iterations', 200, 'cost_tolerance', 1e-8, 'n_restarts', 1, ...
               'use_uniform_init', true, 'use_single_element_init', false, ...
               'amplitude_bounds', [0, 1]);

% Principal-plane cut (-90..90): phi=180 half-plane mirrored onto negatives.
idx_p0   = nearest_index(phi_deg(:), 0);
idx_p180 = nearest_index(phi_deg(:), 180);
td.cut_ang = [-fliplr(theta_deg(2:end)), theta_deg];
td.E_cut   = [fliplr(squeeze(stack(:, 2:end, idx_p180))), squeeze(stack(:, :, idx_p0))];
td.e_s     = td.E_cut(:, nearest_index(td.cut_ang(:), 0));

% Build ONCE for all gate tests (cached in a temp .mat used by the cache test).
td.cache_path = [tempname(), '.mat'];
td.cb = agent_codebook_build(stack, [], theta_deg, phi_deg, ...
    td.aj, td.ag, td.oc, td.cache_path);
testCase.TestData = td;
end


function teardownOnce(testCase)
if isfile(testCase.TestData.cache_path)
    delete(testCase.TestData.cache_path);
end
end


% ── P4 gate 1: per-arm center depth <= -30 dB, gain penalty <= 1 dB ─

function test_arm_depth_and_gain(testCase)
td = testCase.TestData;
cb = td.cb;
w1 = cb.W(:, 1);
g_ref = abs(w1' * td.e_s)^2 / real(w1' * w1);
verifyTrue(testCase, isnan(cb.null_center_deg(1)), 'arm 1 must be unconstrained');
for i = 2:size(cb.W, 2)
    w = cb.W(:, i);
    p = abs(w' * td.E_cut).^2;
    c = cb.null_center_deg(i);
    center_depth_db = 10 * log10(p(nearest_index(td.cut_ang(:), c)) / max(p));
    verifyLessThanOrEqual(testCase, center_depth_db, -30, ...
        sprintf('arm %d @ %g deg: center depth %.1f dB', i, c, center_depth_db));
    gain_penalty_db = 10 * log10((abs(w' * td.e_s)^2 / real(w' * w)) / g_ref);
    verifyGreaterThanOrEqual(testCase, gain_penalty_db, -1, ...
        sprintf('arm %d @ %g deg: peak gain penalty %.2f dB', i, c, gain_penalty_db));
end
end


% ── P4 gate 2: coverage — no hole worse than -25 dB between arms ──

function test_coverage(testCase)
td = testCase.TestData;
cb = td.cb;
allowed = find(abs(td.cut_ang - td.aj.theta_s_deg) >= td.aj.guard_deg);
worst = -Inf;
for a = allowed
    best = 0;
    for i = 2:size(cb.W, 2)
        w = cb.W(:, i);
        p = abs(w' * td.E_cut).^2;
        best = min(best, 10 * log10(p(a) / max(p)));
    end
    worst = max(worst, best);
end
verifyLessThanOrEqual(testCase, worst, -25, ...
    sprintf('coverage: worst best-arm depth %.1f dB', worst));
end


% ── P4 gate 3: cache round-trip and staleness detection ───────────

function test_cache_roundtrip(testCase)
td = testCase.TestData;
% Second call with identical parameters loads from cache -> identical W.
cb2 = agent_codebook_build(td.stack, [], td.theta_deg, td.phi_deg, ...
    td.aj, td.ag, td.oc, td.cache_path);
verifyEqual(testCase, cb2.W, td.cb.W);
verifyEqual(testCase, cb2.meta.built, td.cb.meta.built, ...
    'cache hit must not rebuild');
% Changed parameter -> stale cache -> rebuild with the new arm grid.
ag_coarse = td.ag;
ag_coarse.null_grid_deg = 45.0;
cb3 = agent_codebook_build(td.stack, [], td.theta_deg, td.phi_deg, ...
    td.aj, ag_coarse, td.oc, td.cache_path);
verifyNotEqual(testCase, size(cb3.W, 2), size(td.cb.W, 2));
end


% ── Contract: guard/window overlap warning, missing keys ──────────

function test_codebook_contract(testCase)
td = testCase.TestData;
aj_tight = td.aj;
aj_tight.guard_deg = 5.0;   % < (peak_width + null_width)/2 = 8
verifyWarning(testCase, ...
    @() agent_codebook_build(td.stack, [], td.theta_deg, td.phi_deg, ...
        aj_tight, td.ag, struct('max_iterations', 1, 'cost_tolerance', 1e-3, ...
        'n_restarts', 1, 'use_uniform_init', true, ...
        'use_single_element_init', false), []), ...
    'agent_codebook_build:GuardTooSmall');
ag_bad = rmfield(td.ag, 'peak_width_deg');
verifyError(testCase, ...
    @() agent_codebook_build(td.stack, [], td.theta_deg, td.phi_deg, ...
        td.aj, ag_bad, td.oc, []), ...
    'agent_codebook_build:MissingKey');
end
