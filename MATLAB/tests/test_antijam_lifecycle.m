function tests = test_antijam_lifecycle
% TEST_ANTIJAM_LIFECYCLE  Null-lifecycle scenario (S6): jammer silent for
% 1 min, ON at a static angle for 1 min, silent again (power mode 'window',
% per-scenario duration override). Verifies the requested behavior:
%   - the array steers a deep null onto the jammer when it appears, and
%   - discovers its disappearance and returns to the jammer-free "perfect"
%     peak steering (quiescent gain) afterwards.
% Observed reference behavior (16-el ULA toy): LCMV re-forms the quiescent
% beam within ~1 s of turn-off (exponential forgetting flushes the jammer
% from R_hat); the bandit returns to NEAR-optimal steering but keeps
% wandering between arms — with the jammer silent all arms' rewards tie
% within ~0.3 dB, so Thompson sampling has no pull toward arm 1. Its "return
% to peak" is therefore asserted on gain penalty, not on arm identity.
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
    stack(n, :, :) = exp(1i * pi * (n - 1 - (n_el - 1) / 2) * u);
end

td = struct();
td.aj = struct('theta_s_deg', 0.0, 'guard_deg', 30.0, 'jn_ratio_db', 20.0, ...
               'sigma_s_db', 0.0, 'cut_type', 'theta_cut', 'cut_phi_deg', 0.0);
td.ag = struct('null_grid_deg', 5.0, 'null_width_deg', 10.0, ...
               'peak_width_deg', 6.0, 'null_weight', 100.0, ...
               'method', 'thompson', 'discount', 0.97, 'window', 50, ...
               'sigma_tilde_db', 1.0);
td.acfg = struct('forgetting_lambda', 0.98, 'diagonal_loading_db', 10);
oc = struct('max_iterations', 200, 'cost_tolerance', 1e-8, 'n_restarts', 1, ...
            'use_uniform_init', true, 'use_single_element_init', false, ...
            'amplitude_bounds', [0, 1]);
td.cb = agent_codebook_build(stack, [], theta_deg, phi_deg, td.aj, td.ag, oc, []);

idx_p0   = nearest_index(phi_deg(:), 0);
idx_p180 = nearest_index(phi_deg(:), 180);
td.cut_ang = [-fliplr(theta_deg(2:end)), theta_deg];
td.E_cut   = [fliplr(squeeze(stack(:, 2:end, idx_p180))), squeeze(stack(:, :, idx_p0))];
td.e_s     = td.E_cut(:, nearest_index(td.cut_ang(:), 0));
w_ref      = adapt_lcmv(eye(n_el), td.e_s, 0);
td.g_ref   = abs(w_ref' * td.e_s)^2 / real(w_ref' * w_ref);

td.scn_cfg = struct('id', 'S6', 'motion', 'static', 'power', 'window', ...
                    'on_time_s', 60.0, 'off_time_s', 120.0, 'duration_s', 180.0);
td.sc  = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, ...
                'seed', 1234);
td.thr_db = 5.0;
testCase.TestData = td;
end


% ── Scenario shape: window power mode + duration override ─────────

function test_window_scenario_shape(testCase)
td  = testCase.TestData;
scn = sim_scenario(td.scn_cfg, td.aj, td.sc);
verifyEqual(testCase, numel(scn.t_s), 3601, 'duration_s override -> 180 s');
verifyFalse(testCase, scn.jammer_on(1));
verifyFalse(testCase, scn.jammer_on(end));
verifyTrue(testCase, all(scn.jammer_on(scn.t_s >= 60 & scn.t_s < 120)));
verifyTrue(testCase, ~any(scn.jammer_on(scn.t_s < 60 | scn.t_s >= 120)));
verifyEqual(testCase, numel(scn.theta_j_deg), 3601);
verifyEqual(testCase, max(scn.theta_j_deg) - min(scn.theta_j_deg), 0, ...
    'static jammer position');
verifyEqual(testCase, cellfun(@(e) e.type, scn.events(:).', 'UniformOutput', false), ...
    {'turn_on', 'turn_off'});
% Missing window keys raise (no silent defaults).
bad = rmfield(td.scn_cfg, 'off_time_s');
verifyError(testCase, @() sim_scenario(bad, td.aj, td.sc), 'sim_scenario:MissingKey');
end


% ── LCMV tracker: null forms on turn-on, beam restored after off ──

function test_lcmv_null_lifecycle(testCase)
td  = testCase.TestData;
scn = sim_scenario(td.scn_cfg, td.aj, td.sc);
st  = sim_engine_init(td.E_cut, [], td.cut_ang, scn, td.aj, td.sc, 'C');
trk = adapt_tracking_init(td.acfg, st.e_s, size(td.E_cut, 1));
T = numel(scn.t_s);
sinr = zeros(1, T); gpen = zeros(1, T); ndepth = NaN(1, T);
jidx = nearest_index(td.cut_ang(:), scn.theta_j_deg(1));
w = trk.w;
for k = 1:T
    gpen(k) = 10 * log10((abs(w' * td.e_s)^2 / real(w' * w)) / td.g_ref);
    p = abs(w' * td.E_cut).^2;
    ndepth(k) = 10 * log10(p(jidx) / max(p));
    [obs, st] = sim_engine_step(st, w);
    sinr(k) = obs.sinr_db;
    [w, trk] = adapt_tracking_update(trk, obs);
end
k_on  = find(scn.t_s >= 60, 1);
k_off = find(scn.t_s >= 120, 1);
on_ss = (k_on + 100):(k_off - 1);          % ON steady state (5 s settle)
post  = (k_off + 100):T;                   % > 5 s after turn-off

% Null steered onto the jammer while it transmits.
verifyLessThanOrEqual(testCase, find(sinr(k_on:end) >= td.thr_db, 1) - 1, 25, ...
    'recovery after turn-on');
verifyGreaterThanOrEqual(testCase, mean(sinr(on_ss) >= td.thr_db), 0.99);
verifyLessThanOrEqual(testCase, median(ndepth(on_ss)), -30, ...
    sprintf('ON null depth median %.1f dB', median(ndepth(on_ss))));

% Jammer disappearance discovered: quiescent beam restored (forgetting
% flushes the jammer from R_hat within ~1 s; assert with margin at 5 s).
verifyLessThanOrEqual(testCase, max(abs(gpen(post))), 0.5, ...
    sprintf('post-off gain penalty max %.2f dB', max(abs(gpen(post)))));
end


% ── Bandit: instant nulling, near-optimal steering after off ──────

function test_bandit_null_lifecycle(testCase)
td  = testCase.TestData;
scn = sim_scenario(td.scn_cfg, td.aj, td.sc);
st  = sim_engine_init(td.E_cut, [], td.cut_ang, scn, td.aj, td.sc, 'S');
bd  = agent_bandit_init(td.ag, td.cb, 42);
T = numel(scn.t_s);
sinr = zeros(1, T); gpen = zeros(1, T);
w = bd.w;
for k = 1:T
    gpen(k) = 10 * log10((abs(w' * td.e_s)^2 / real(w' * w)) / td.g_ref);
    [obs, st] = sim_engine_step(st, w);
    sinr(k) = obs.sinr_db;
    [w, bd] = agent_bandit_update(bd, obs);
end
k_on  = find(scn.t_s >= 60, 1);
k_off = find(scn.t_s >= 120, 1);
on_ss = (k_on + 100):(k_off - 1);
post  = (k_off + 100):T;

verifyLessThanOrEqual(testCase, find(sinr(k_on:end) >= td.thr_db, 1) - 1, 30, ...
    'recovery after turn-on');
verifyGreaterThanOrEqual(testCase, mean(sinr(on_ss) >= td.thr_db), 0.95);
% Near-optimal steering after the jammer disappears: all arms tie within
% ~0.3 dB, so assert gain penalty, not arm identity (see file header).
verifyLessThanOrEqual(testCase, mean(abs(gpen(post))), 1.0, ...
    sprintf('post-off mean gain penalty %.2f dB', mean(abs(gpen(post)))));
end
