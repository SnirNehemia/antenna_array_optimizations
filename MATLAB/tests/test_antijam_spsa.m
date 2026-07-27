function tests = test_antijam_spsa
% TEST_ANTIJAM_SPSA  P3 gate for the Mode S SPSA baseline
% (antijam_milestone_plan.md Section 4): static jammer, reach within 3 dB of
% the oracle SINR in <= 150 SINR probes, median over 50 Monte Carlo seeds.
% Tuning per the P3 sweep: a = 2, c = 0.2, A = 15, step_max = 0.3 — the
% stability constant and step clamp are load-bearing (see adapt_spsa_init).
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));

testCase.TestData.aj = struct('theta_s_deg', 0.0, 'phi_s_deg', 0.0, ...
    'guard_deg', 10.0, 'jn_ratio_db', 20.0, 'sigma_s_db', 0.0);
testCase.TestData.spsa_cfg = struct('a', 2.0, 'c', 0.2, 'alpha', 0.602, ...
    'gamma', 0.101, 'A', 15, 'step_max', 0.3);
testCase.TestData.theta_deg = linspace(0, 180, 361);
testCase.TestData.phi_deg   = 0;
testCase.TestData.E   = exp(1i * pi * (0:7).' * sind(testCase.TestData.theta_deg));
testCase.TestData.stack = reshape(testCase.TestData.E, 8, 361, 1);
end


% ── P3 gate: within 3 dB of oracle in <= 150 probes (median, 50 MC) ─

function test_static_convergence_monte_carlo(testCase)
td = testCase.TestData;
n_seeds  = 50;
n_probes = 150;
hit = inf(1, n_seeds);
for s = 1:n_seeds
    sc  = struct('dt_s', 0.05, 'duration_s', 60.0, 'snapshots_per_step', 16, ...
                 'seed', 1000 + s);
    scn = sim_scenario(struct('id', 'S1', 'motion', 'static', 'power', 'constant'), ...
                       td.aj, sc);

    % Oracle SINR for this seed's jammer angle (static -> constant over time).
    st = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, 'S');
    [~, st1] = sim_engine_step(st, ones(8, 1));
    R = sim_analytic_covariance(st1);
    oracle_db = 10 * log10(st1.sigma_s_sq * real(st1.e_s' * (R \ st1.e_s)));

    % SPSA closed loop from uniform weights.
    st   = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, 'S');
    spsa = adapt_spsa_init(td.spsa_cfg, ones(8, 1), 2000 + s);
    w = spsa.w;
    for k = 1:n_probes
        [obs, st] = sim_engine_step(st, w);
        if obs.sinr_db >= oracle_db - 3
            hit(s) = k;
            break;
        end
        [w, spsa] = adapt_spsa_update(spsa, obs);
    end
end
verifyLessThanOrEqual(testCase, median(hit), 150, ...
    sprintf('median probes to oracle-3dB = %g', median(hit)));
% Robustness margin beyond the formal gate: most seeds must actually get there.
verifyGreaterThanOrEqual(testCase, nnz(isfinite(hit)), 45, ...
    sprintf('%d/%d seeds converged within %d probes', ...
    nnz(isfinite(hit)), n_seeds, n_probes));
end


% ── Contract checks ───────────────────────────────────────────────

function test_spsa_contract(testCase)
td = testCase.TestData;
% Missing required config key raises (no silent defaults).
bad_cfg = rmfield(td.spsa_cfg, 'step_max');
verifyError(testCase, @() adapt_spsa_init(bad_cfg, ones(8, 1), 1), ...
    'adapt_spsa_init:MissingKey');

% Probe pairing: first two returned weight vectors are x +/- c*delta.
spsa = adapt_spsa_init(td.spsa_cfg, ones(8, 1), 1);
w_plus = spsa.w;
[w_minus, spsa] = adapt_spsa_update(spsa, struct('sinr_db', 0, 'snapshots', []));
x0 = weights_to_x(ones(8, 1));
verifyEqual(testCase, weights_to_x(w_plus) + weights_to_x(w_minus), 2 * x0, ...
    'AbsTol', 1e-12, 'probe pair must be symmetric about x_est');

% Deterministic under seed.
spsa_a = adapt_spsa_init(td.spsa_cfg, ones(8, 1), 42);
spsa_b = adapt_spsa_init(td.spsa_cfg, ones(8, 1), 42);
verifyEqual(testCase, spsa_a.w, spsa_b.w);
end
