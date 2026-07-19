function tests = test_antijam_bandit
% TEST_ANTIJAM_BANDIT  P5 gates for the non-stationary bandit agent
% (antijam_milestone_plan.md Section 4) on the 16-el ULA + P4 codebook:
%   1. Static jammer: within 30 probes the played arm covers theta_j OR is
%      within 1 dB of the best arm's SINR, >= 90/100 Monte Carlo runs. (The
%      1 dB clause refines the plan wording: arms whose natural sidelobe
%      nulls coincide with theta_j earn rewards identical to the designated
%      arm, so pure arm identification is ill-posed — SINR-optimality is the
%      operational criterion.)
%   2. Drifting jammer (S2): SINR availability >= 90%.
%   3. Recovery after the S3 10-deg jump: bandit beats SPSA by >= 2x (median).
% Tuning: SIGMA_TILDE = 1 dB (see agent_bandit_init), discount 0.97.
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
               'method', 'thompson', 'discount', 0.97, 'window', 50);
td.spsa_cfg = struct('a', 2.0, 'c', 0.2, 'alpha', 0.602, 'gamma', 0.101, ...
                     'A', 15, 'step_max', 0.3);
oc = struct('max_iterations', 200, 'cost_tolerance', 1e-8, 'n_restarts', 1, ...
            'use_uniform_init', true, 'use_single_element_init', false, ...
            'amplitude_bounds', [0, 1]);

td.cb = agent_codebook_build(stack, [], theta_deg, phi_deg, td.aj, td.ag, oc, []);

idx_p0   = nearest_index(phi_deg(:), 0);
idx_p180 = nearest_index(phi_deg(:), 180);
td.cut_ang = [-fliplr(theta_deg(2:end)), theta_deg];
td.E_cut   = [fliplr(squeeze(stack(:, 2:end, idx_p180))), squeeze(stack(:, :, idx_p0))];
td.thr_db  = 5.0;
testCase.TestData = td;
end


% ────────────────────────── HELPERS ───────────────────────────────

function sinr_arm = per_arm_sinr(td, theta_j_deg)
% True SINR of every codebook arm against a jammer at theta_j (J/N 20 dB).
aidx = nearest_index(td.cut_ang(:), theta_j_deg);
sidx = nearest_index(td.cut_ang(:), td.aj.theta_s_deg);
n_arms = size(td.cb.W, 2);
sinr_arm = zeros(1, n_arms);
for i = 1:n_arms
    wi = td.cb.W(:, i);
    gs = abs(wi' * td.E_cut(:, sidx))^2;
    gj = abs(wi' * td.E_cut(:, aidx))^2;
    sinr_arm(i) = 10 * log10(gs / (100 * gj + real(wi' * wi)));
end
end


function [sinr, bd] = bandit_loop(td, scn, seed, n_probes)
sc = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, 'seed', seed);
st = sim_engine_init(td.E_cut, [], td.cut_ang, scn, td.aj, sc, 'S');
bd = agent_bandit_init(td.ag, td.cb, seed + 2000);
w = bd.w;
sinr = zeros(1, n_probes);
for k = 1:n_probes
    [obs, st] = sim_engine_step(st, w);
    sinr(k) = obs.sinr_db;
    [w, bd] = agent_bandit_update(bd, obs);
end
end


% ── P5 gate 1: static identification within 30 probes, 90/100 ─────

function test_static_identification(testCase)
td = testCase.TestData;
ok = 0;
for s = 1:100
    sc  = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, ...
                 'seed', 1000 + s);
    scn = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), ...
                       td.aj, sc);
    [~, bd] = bandit_loop(td, scn, 1000 + s, 30);
    c = td.cb.null_center_deg(bd.arm_index);
    covered = ~isnan(c) && abs(c - scn.theta_j_deg(1)) <= td.ag.null_width_deg;
    sinr_arm = per_arm_sinr(td, scn.theta_j_deg(1));
    near_opt = (max(sinr_arm) - sinr_arm(bd.arm_index)) <= 1;
    ok = ok + double(covered || near_opt);
end
verifyGreaterThanOrEqual(testCase, ok, 90, ...
    sprintf('static identification: %d/100 within 30 probes', ok));
end


% ── P5 gate 2: drifting jammer availability >= 90% ────────────────

function test_drift_availability(testCase)
td = testCase.TestData;
avail = zeros(1, 10);
for s = 1:10
    sc  = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, ...
                 'seed', 1200 + s);
    scn = sim_scenario(struct('id', 'S2', 'motion', 'drift', ...
                'drift_deg_per_s', 2.0, 'power', 'constant'), td.aj, sc);
    sinr = bandit_loop(td, scn, 1200 + s, numel(scn.t_s));
    avail(s) = mean(sinr >= td.thr_db);
end
verifyGreaterThanOrEqual(testCase, median(avail), 0.90, ...
    sprintf('S2 availability median %.3f (min %.3f)', median(avail), min(avail)));
end


% ── P5 gate 3: recovery after jump beats SPSA by >= 2x (median) ───

function test_recovery_beats_spsa(testCase)
td = testCase.TestData;
n_seeds = 10;
rec_b = zeros(1, n_seeds);
rec_s = zeros(1, n_seeds);
for s = 1:n_seeds
    sc  = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, ...
                 'seed', 1500 + s);
    scn = sim_scenario(struct('id', 'S3', 'motion', 'drift', 'drift_deg_per_s', 2.0, ...
                'jump_deg', 10.0, 'jump_time_s', 30.0, 'power', 'constant'), td.aj, sc);
    k_jump = find(scn.t_s >= 30.0, 1);
    T = numel(scn.t_s);

    sinr = bandit_loop(td, scn, 1500 + s, T);
    r = find(sinr(k_jump:end) >= td.thr_db, 1) - 1;
    if isempty(r), r = Inf; end
    rec_b(s) = r;

    st = sim_engine_init(td.E_cut, [], td.cut_ang, scn, td.aj, sc, 'S');
    sp = adapt_spsa_init(td.spsa_cfg, ones(size(td.E_cut, 1), 1), 1500 + s);
    w = sp.w;
    sinr = zeros(1, T);
    for k = 1:T
        [obs, st] = sim_engine_step(st, w);
        sinr(k) = obs.sinr_db;
        [w, sp] = adapt_spsa_update(sp, obs);
    end
    r = find(sinr(k_jump:end) >= td.thr_db, 1) - 1;
    if isempty(r), r = Inf; end
    rec_s(s) = r;
end
% max(.,1) keeps the ratio well-defined when the bandit recovers instantly.
verifyGreaterThanOrEqual(testCase, median(rec_s), 2 * max(median(rec_b), 1), ...
    sprintf('recovery medians: bandit %g, SPSA %g', median(rec_b), median(rec_s)));
end


% ── swucb alternative: smoke (converges on a static jammer) ───────

function test_swucb_smoke(testCase)
td = testCase.TestData;
td.ag.method = 'swucb';
sc  = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, 'seed', 1042);
scn = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), ...
                   td.aj, sc);
sinr = bandit_loop(td, scn, 1042, 400);
verifyGreaterThan(testCase, mean(sinr(201:end) >= td.thr_db), 0.7, ...
    'swucb should mostly exploit good arms after burn-in');
end


% ── Contract checks ───────────────────────────────────────────────

function test_bandit_contract(testCase)
td = testCase.TestData;
ag_bad = td.ag;
ag_bad.method = 'greedy';
verifyError(testCase, @() agent_bandit_init(ag_bad, td.cb, 1), ...
    'agent_bandit_init:BadMethod');
verifyError(testCase, @() agent_bandit_init(rmfield(td.ag, 'discount'), td.cb, 1), ...
    'agent_bandit_init:MissingKey');
% Deterministic under seed.
b1 = agent_bandit_init(td.ag, td.cb, 7);
b2 = agent_bandit_init(td.ag, td.cb, 7);
obs = struct('sinr_db', 3.0, 'snapshots', []);
[w1, ~] = agent_bandit_update(b1, obs);
[w2, ~] = agent_bandit_update(b2, obs);
verifyEqual(testCase, w1, w2);
end
