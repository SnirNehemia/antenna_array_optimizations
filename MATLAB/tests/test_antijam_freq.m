function tests = test_antijam_freq
% TEST_ANTIJAM_FREQ  P10 gates for jammer carrier estimation + RF notch modeling
% (antijam_milestone_plan.md Section 4, P10).
%
% The phase adds a TEMPORAL axis to a simulator that had none: pre-P10 the
% snapshot columns were i.i.d. Gaussian draws, spectrally white by construction,
% so there was no frequency to estimate. The two gates that actually protect the
% existing milestone are therefore the compatibility ones:
%
%   G1 (byte-identity) — with sim.fs_hz absent, the snapshot generator must be
%      bit-for-bit what it was before, floating-point rounding included. Checked
%      against an independent inline reimplementation of the pre-P10 formula,
%      not against a recorded fixture, so it cannot rot into a tautology.
%   G2 (covariance invariance) — with the waveform layer ON, the block-averaged
%      SPATIAL covariance must still converge to sim_analytic_covariance. Every
%      Mode C algorithm consumes only (X*X')/K, so if this holds none of them can
%      tell a CW tone from Gaussian noise. G7 confirms that end-to-end on the P2
%      tracker itself.
%
% The remaining gates cover the estimator (G3), the notch model's closed forms
% (G4, G5), the closed loop through kpi_evaluate (G6), no regression on the P2
% tracker (G7), config/contract validation (G8), the two cases where the notch
% actually earns its place (G9), and hold-through-OFF (G10).
%
% Same toy 8-element half-wavelength ULA as test_antijam_tracking.m, but with
% K = 128 snapshots/step (the value config.yaml actually uses, and the one the
% periodogram resolution depends on) rather than that suite's K = 16.
%
% CAVEAT ON THIS FIXTURE, worth knowing before reading G6's "the notch adds
% nothing" bound as a general claim: an 8-element single-polarization ULA facing
% ONE jammer has a large surplus of degrees of freedom, so its spatial null is
% essentially perfect and leaves the notch nothing to remove. The real 6-element
% dual-polarization array is DOF-limited (rank-2 desired + rank-2 jammer out of
% 6), its null is much shallower, and there the notch is worth 9-11 dB even
% off-beam — see results/freq_notch_demo/ and the [P10] plan section. The gates
% here bound the mechanism; the demo measures the payoff.
%
% Thresholds are calibrated from measured sweeps, not guessed — see the inline
% notes at each gate.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
addpath(fullfile(fileparts(here), 'antijam_utils'));

theta_deg = linspace(0, 180, 361);
phi_deg   = 0;
E     = exp(1i * pi * (0:7).' * sind(theta_deg));   % half-wavelength 8-el ULA
stack = reshape(E, 8, 361, 1);

testCase.TestData.theta_deg = theta_deg;
testCase.TestData.phi_deg   = phi_deg;
testCase.TestData.E         = E;
testCase.TestData.stack     = stack;
testCase.TestData.aj = struct('theta_s_deg', 0.0, 'phi_s_deg', 0.0, 'guard_deg', 10.0, ...
                              'jn_ratio_db', 20.0, 'sigma_s_db', 0.0, 'f_center_hz', 2.4e9);
% Base sim config WITHOUT fs_hz — the pre-P10 path. Gates add fs_hz explicitly.
testCase.TestData.sc = struct('dt_s', 0.05, 'duration_s', 10.0, ...
                              'snapshots_per_step', 128, 'seed', 1234);
testCase.TestData.fs_hz     = 24e6;
testCase.TestData.freq_cfg  = struct('nfft_factor', 4, 'presence_snr_db', 15.0);
testCase.TestData.notch_cfg = struct('mode', 'adaptive', 'depth_db', 35.0, ...
                                     'bw_hz', 2.0e5, 'adaptation_tap', 'pre');
end


% ────────────────────────── HELPERS ───────────────────────────────

function scn = make_scn(td, varargin)
% Static jammer at 40 deg, constant power; extra scenario keys passed through.
cfg = struct('id', 'F1', 'motion', 'static', 'power', 'constant', ...
             'theta_j_deg', 40.0, 'phi_j_deg', 0.0);
for i = 1:2:numel(varargin)
    cfg.(varargin{i}) = varargin{i + 1};
end
scn = sim_scenario(cfg, td.aj, td.sc);
end


function sc = with_fs(td)
sc = td.sc;
sc.fs_hz = td.fs_hz;
end


function gap = tracker_gap(td, sc, scn)
% Mean steady-state oracle gap of the P2 covariance tracker (second half of the
% run), the same honest closed loop test_antijam_tracking.m uses.
acfg = struct('forgetting_lambda', 0.90, 'diagonal_loading_db', 10);
st   = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, 'C');
trk  = adapt_tracking_init(acfg, st.e_s, size(td.E, 1));
T    = numel(scn.t_s);
g    = zeros(1, T);
w    = trk.w;
for k = 1:T
    [obs, st] = sim_engine_step(st, w);
    R = sim_analytic_covariance(st);
    g(k) = 10 * log10(st.sigma_s_sq * real(st.e_s' * (R \ st.e_s))) - obs.sinr_db;
    [w, trk] = adapt_tracking_update(trk, obs);
end
gap = mean(g(floor(T / 2):end));
end


% ────────────────────────── GATES ─────────────────────────────────

function test_g1_byte_identical_without_fs(testCase)
% G1: sim.fs_hz absent -> snapshots and SINR bit-for-bit identical to the
% pre-P10 generator. Reference is reimplemented inline from the old formula
% (signal -> jammer -> noise, i.i.d. cgauss, that exact accumulation order).
td  = testCase.TestData;
scn = make_scn(td);
w   = ones(8, 1) / 8;
st  = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, td.sc, 'C');

ref = RandStream('mt19937ar', 'Seed', td.sc.seed);
cg  = @(s, m, n, v) sqrt(v / 2) * (randn(s, m, n) + 1i * randn(s, m, n));
it  = nearest_index_2d(td.theta_deg, td.phi_deg, scn.theta_j_deg(1), scn.phi_j_deg(1));
e_j = td.E(:, it);
sj  = 10^(td.aj.jn_ratio_db / 10);
K   = td.sc.snapshots_per_step;

for k = 1:5
    [obs, st] = sim_engine_step(st, w);
    xr = st.e_s * cg(ref, 1, K, 10^(td.aj.sigma_s_db / 10));
    xr = xr + e_j * cg(ref, 1, K, sj);
    xr = xr + cg(ref, 8, K, 1.0);
    verifyEqual(testCase, obs.snapshots, xr, ...
        'G1: snapshots differ from the pre-P10 generator (must be bit-identical).');
    verifyEmpty(testCase, obs.fs_hz, 'G1: obs.fs_hz must stay empty with the waveform layer off.');
end
verifyFalse(testCase, st.use_waveform, 'G1: use_waveform must be false without sim.fs_hz.');
end


function test_g2_spatial_covariance_invariance(testCase)
% G2: with the CW tone on, the block-averaged spatial covariance still converges
% to the analytic R. Same 5% Frobenius threshold as the P1 convergence gate.
% Measured on this fixture: ~1% at M*K = 6400 (= 800 x N_el).
td  = testCase.TestData;
sc  = with_fs(td);
scn = make_scn(td, 'freq', 'constant', 'f_j_offset_hz', 3.1e6);
st  = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, 'C');
w   = ones(8, 1) / 8;
M   = 50;
R   = zeros(8);
for k = 1:M
    [obs, st] = sim_engine_step(st, w);
    R = R + (obs.snapshots * obs.snapshots') / size(obs.snapshots, 2);
end
R  = R / M;
Ra = sim_analytic_covariance(st);
rel = norm(R - Ra, 'fro') / norm(Ra, 'fro');
fprintf('  G2 covariance rel Frobenius error = %.2f%%\n', 100 * rel);
verifyLessThan(testCase, rel, 0.05, ...
    'G2: tone-driven snapshots must still give the analytic spatial covariance.');
verifyTrue(testCase, st.use_waveform);
verifyEqual(testCase, obs.fs_hz, td.fs_hz);
end


function test_g3_estimator_accuracy(testCase)
% G3: single-block carrier accuracy over 100 random carriers. Thresholds
% calibrated from the measured sweep (see the fprintf output): what limits
% accuracy is the in-band ratio sigma_j^2/(sigma_s^2 + sigma_n^2) — the WIDEBAND
% DESIRED SIGNAL, not the noise floor — so the bound is checked at both the
% quiet (sigma_s_db = 0) and loud (sigma_s_db = 20) operating points.
td = testCase.TestData;
sc = with_fs(td);
for sigma_s_db = [0, 20]
    aj = td.aj;
    aj.sigma_s_db = sigma_s_db;
    nT  = 100;
    err = NaN(1, nT);
    rng(7);
    for t = 1:nT
        f_true = (rand - 0.5) * 0.8 * td.fs_hz;
        s_cfg  = sc;
        s_cfg.seed = 1000 + t;
        cfg = struct('id', 'G3', 'motion', 'static', 'power', 'constant', ...
                     'theta_j_deg', 40.0, 'phi_j_deg', 0.0, ...
                     'freq', 'constant', 'f_j_offset_hz', f_true);
        scn = sim_scenario(cfg, aj, s_cfg);
        st  = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, aj, s_cfg, 'C');
        [obs, ~] = sim_engine_step(st, ones(8, 1) / 8);
        fq = adapt_freq_estimate(obs.snapshots, obs.fs_hz, td.freq_cfg);
        verifyTrue(testCase, fq.present, ...
            sprintf('G3: missed detection at sigma_s_db = %d (J/N = 20 dB).', sigma_s_db));
        err(t) = fq.f_hat_hz - f_true;
    end
    rmse = sqrt(mean(err.^2));
    fprintf('  G3 sigma_s_db=%2d: RMSE %6.0f Hz, max |e| %6.0f Hz (bin = %.0f Hz)\n', ...
        sigma_s_db, rmse, max(abs(err)), td.fs_hz / (td.freq_cfg.nfft_factor * 128));
    % 5 kHz is ~8x the measured RMSE at the harder operating point and still
    % only 2.5% of a 200 kHz notch — a real bound, not a rubber stamp.
    verifyLessThan(testCase, rmse, 5e3, ...
        sprintf('G3: carrier RMSE too high at sigma_s_db = %d.', sigma_s_db));
end
end


function test_g4_notch_closed_forms(testCase)
% G4: the notch response matches its own definition, and mode 'off' is inert.
td = testCase.TestData;
n  = td.notch_cfg;
d  = 10^(-n.depth_db / 10);

[h_centre, loss] = sim_notch_response(n, 0, 0, td.fs_hz);
verifyEqual(testCase, h_centre, d, 'RelTol', 1e-12, ...
    'G4: |H|^2 at the notch centre must equal the configured depth.');

h_edge = sim_notch_response(n, 0, n.bw_hz / 2, td.fs_hz);
verifyEqual(testCase, h_edge, (1 + d) / 2, 'RelTol', 1e-12, ...
    'G4: bw_hz must be the true 3-dB width.');

loss_ref = pi * (1 - d) * n.bw_hz / (2 * td.fs_hz);
fprintf('  G4 loss_frac %.5f vs closed form %.5f (%.3f dB signal cost)\n', ...
    loss, loss_ref, -10 * log10(1 - loss));
verifyEqual(testCase, loss, loss_ref, 'RelTol', 0.01, ...
    'G4: wideband insertion loss must match pi*(1-d)*B/(2*fs) for a narrow notch.');

% mode 'off' must reproduce the un-notched SINR exactly.
sc  = with_fs(td);
scn = make_scn(td, 'freq', 'constant', 'f_j_offset_hz', 3.1e6);
n_off = n; n_off.mode = 'off';
st  = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, 'C', n_off);
[obs, st] = sim_engine_step(st, ones(8, 1) / 8, 3.1e6);
verifyEqual(testCase, obs.sinr_db, st.last.sinr_no_notch_db, 'AbsTol', 0, ...
    'G4: notch.mode ''off'' must leave SINR untouched.');
end


function test_g5_misestimation_sensitivity(testCase)
% G5: the SINR benefit the engine delivers is exactly the notch response
% evaluated at the estimation error — this is what makes G3's accuracy KPI
% mean something — and it degrades monotonically as the estimate drifts off.
td    = testCase.TestData;
sc    = with_fs(td);
f_true = 3.1e6;
scn   = make_scn(td, 'freq', 'constant', 'f_j_offset_hz', f_true);
offsets = [0, 1e3, 1e4, 5e4, 1e5, 5e5, 5e6];
gain_db = zeros(size(offsets));
for i = 1:numel(offsets)
    st = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, 'C', td.notch_cfg);
    [obs, st] = sim_engine_step(st, ones(8, 1) / 8, f_true + offsets(i));
    gain_db(i) = obs.sinr_db - st.last.sinr_no_notch_db;
    h2_expect  = sim_notch_response(td.notch_cfg, f_true + offsets(i), f_true, td.fs_hz);
    verifyEqual(testCase, st.last.notch_rejection_db, -10 * log10(h2_expect), 'RelTol', 1e-12, ...
        'G5: engine rejection must equal sim_notch_response at the estimation error.');
end
fprintf('  G5 SINR gain vs |f_hat - f_j|: %s dB  (at offsets %s kHz)\n', ...
    num2str(gain_db, '%.1f '), num2str(offsets / 1e3, '%g '));
verifyTrue(testCase, all(diff(gain_db) < 1e-9), ...
    'G5: notch benefit must decrease monotonically with carrier-estimation error.');
verifyLessThan(testCase, gain_db(end), 0, ...
    'G5: a notch placed far from the jammer must be a net loss (insertion loss only).');

% The delivered gain must equal the closed form exactly — a far stronger claim
% than any magic threshold, and it pins down the CEILING: a notch can only
% remove the jammer's contribution, so its benefit saturates at the
% jammer-to-noise ratio AT THE BEAMFORMER OUTPUT. With this quiescent uniform
% beam the array's own sidelobe already suppresses a 40 deg jammer by ~17 dB,
% which is why a 35 dB notch buys only ~12 dB here. See G9: the gain grows to
% ~28 dB when the jammer sits where the pattern cannot suppress it at all.
[~, loss] = sim_notch_response(td.notch_cfg, f_true, f_true, td.fs_hz);
h2_0  = sim_notch_response(td.notch_cfg, f_true, f_true, td.fs_hz);
w     = ones(8, 1) / 8;
st0   = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, 'C', td.notch_cfg);
it    = nearest_index_2d(td.theta_deg, td.phi_deg, scn.theta_j_deg(1), scn.phi_j_deg(1));
p_sig = st0.sigma_s_sq * abs(w' * st0.e_s)^2;
p_jam = 10^(td.aj.jn_ratio_db / 10) * abs(w' * td.E(:, it))^2;
p_nse = real(w' * w);
expect = 10 * log10(p_sig * (1 - loss) / (p_jam * h2_0 + p_nse * (1 - loss))) ...
       - 10 * log10(p_sig / (p_jam + p_nse));
verifyEqual(testCase, gain_db(1), expect, 'RelTol', 1e-10, ...
    'G5: delivered SINR gain must match the closed-form notch model.');
ceiling_db = 10 * log10(1 + p_jam / p_nse);
verifyLessThan(testCase, gain_db(1), ceiling_db, ...
    'G5: notch gain cannot exceed the jammer-to-noise ratio at the array output.');
end


function test_g6_closed_loop_tracker_and_notch(testCase)
% G6: end-to-end through closed_loop_run + kpi_evaluate, on a HOPPING carrier —
% the tracker must re-lock after each hop rather than smoothing slowly across
% the gap (the HOP branch of adapt_freq_update).
%
% This gate also pins the ceiling on what a notch can do: it removes only what
% the jammer still CONTRIBUTES, so its benefit is bounded by the jammer-to-noise
% ratio at the beamformer output. On this DOF-rich toy fixture the LCMV null has
% already driven that to nothing, so the notch measures ~0 dB here — see the
% fixture caveat in the file header before generalizing that. G9 covers the two
% cases where the bound is loose and the notch pays: the main-beam jammer the
% spatial nuller declares out of scope, and the turn-on transient before the
% covariance has adapted.
td  = testCase.TestData;
sc  = with_fs(td);
cfg = struct('adapt', struct('forgetting_lambda', 0.90, 'diagonal_loading_db', 10, ...
                             'freq', struct('nfft_factor', 4, 'presence_snr_db', 15.0, ...
                                            'smoothing_lambda', 0.7)), ...
             'notch', td.notch_cfg);
scn_cfg = struct('id', 'G6', 'motion', 'static', 'power', 'constant', ...
                 'theta_j_deg', 40.0, 'phi_j_deg', 0.0, 'freq', 'hop', ...
                 'f_j_hop_offsets_hz', [3.1e6, -7.45e6], 'f_j_hop_period_s', 4.0);
scn = sim_scenario(scn_cfg, td.aj, sc);
log = closed_loop_run('lcmv', td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, cfg, []);
log.oracle_sinr_db = zeros(1, numel(scn.t_s));   % unused by the P10 KPIs
kpi = kpi_evaluate(log, scn, struct('sinr_min_db', 10.0));

n_hops = sum(strcmp(cellfun(@(e) e.type, scn.events, 'UniformOutput', false), 'freq_hop'));
fprintf('  G6 freq RMSE %.0f Hz | %d hops re-locked | notch gain over a converged null %+.2f dB\n', ...
    kpi.freq_rmse_hz, n_hops, kpi.notch_gain_on_mean_db);
verifyTrue(testCase, all(log.freq_present), 'G6: an always-on jammer must be detected every step.');
verifyGreaterThan(testCase, n_hops, 0, 'G6: the fixture must actually hop.');
verifyLessThan(testCase, kpi.freq_rmse_hz, 5e3, ...
    'G6: closed-loop carrier RMSE too high — the tracker is not re-locking after hops.');
verifyEqual(testCase, numel(log.notch_f_hz), numel(scn.t_s));

% The honest bound, both ways: over a converged spatial null the notch neither
% helps materially nor hurts. If this ever grows, the null has stopped working.
verifyLessThan(testCase, abs(kpi.notch_gain_on_mean_db), 0.5, ...
    'G6: over a converged LCMV null there is no jammer left for the notch to remove.');
verifyGreaterThan(testCase, kpi.notch_gain_on_mean_db, -0.2, ...
    'G6: the notch must never be a net loss while the jammer is transmitting.');
end


function test_g9_notch_covers_the_nulls_blind_spots(testCase)
% G9: where the RF notch actually earns its place. Two cases the spatial null
% cannot handle, measured against the notch-off baseline on the same scenario:
%
%   (a) MAIN-BEAM JAMMER (theta_j = theta_s). Nulling boresight would destroy
%       the desired signal, so LCMV correctly refuses — this is exactly the case
%       antijam.guard_deg exists to exclude ("jammer inside the main beam is out
%       of scope for a single-jammer milestone", plan Section 5). The notch is
%       orthogonal to angle and rescues it outright.
%   (b) TURN-ON TRANSIENT. Immediately after the jammer appears, the covariance
%       estimate has not adapted and the null is not yet formed; the notch is
%       already on the right carrier because adapt_freq_update HELD it through
%       the OFF gap.
td = testCase.TestData;
sc = with_fs(td);
base = struct('adapt', struct('forgetting_lambda', 0.90, 'diagonal_loading_db', 10, ...
                  'freq', struct('nfft_factor', 4, 'presence_snr_db', 15.0, ...
                                 'smoothing_lambda', 0.7)));

% (a) Main-beam jammer: guard_deg = 0 so the scenario generator will place it there.
aj = td.aj;
aj.guard_deg = 0.0;
scn_cfg = struct('id', 'G9a', 'motion', 'static', 'power', 'constant', ...
                 'theta_j_deg', 0.0, 'phi_j_deg', 0.0, ...
                 'freq', 'constant', 'f_j_offset_hz', 3.1e6);
scn  = sim_scenario(scn_cfg, aj, sc);
mean_db = zeros(1, 2);
modes = {'off', 'adaptive'};
for i = 1:2
    cfg = base;
    cfg.notch = td.notch_cfg;
    cfg.notch.mode = modes{i};
    log = closed_loop_run('lcmv', td.stack, [], td.theta_deg, td.phi_deg, scn, aj, sc, cfg, []);
    T = numel(scn.t_s);
    mean_db(i) = mean(log.sinr_db(floor(T / 2):end));
end
fprintf('  G9a main-beam jammer: notch off %+.2f dB -> adaptive %+.2f dB (rescue %.1f dB)\n', ...
    mean_db(1), mean_db(2), mean_db(2) - mean_db(1));
verifyGreaterThan(testCase, mean_db(2) - mean_db(1), 20, ...
    'G9a: the notch must rescue the main-beam jammer the spatial null cannot touch.');

% (b) Turn-on transient, jammer well outside the main beam so the null does work.
scn_cfg = struct('id', 'G9b', 'motion', 'static', 'power', 'onoff', ...
                 'duty_cycle', 0.5, 'toggle_period_s', 4.0, ...
                 'theta_j_deg', 40.0, 'phi_j_deg', 0.0, ...
                 'freq', 'constant', 'f_j_offset_hz', 3.1e6);
scn  = sim_scenario(scn_cfg, td.aj, sc);
first_step = zeros(1, 2);
for i = 1:2
    cfg = base;
    cfg.notch = td.notch_cfg;
    cfg.notch.mode = modes{i};
    log = closed_loop_run('lcmv', td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, cfg, []);
    on  = scn.jammer_on;
    ton = find(on & [false, ~on(1:end - 1)], 1);
    first_step(i) = log.sinr_db(ton);
end
fprintf('  G9b first step after turn-on: notch off %+.2f dB -> adaptive %+.2f dB (%.1f dB of lag removed)\n', ...
    first_step(1), first_step(2), first_step(2) - first_step(1));
verifyGreaterThan(testCase, first_step(2) - first_step(1), 3, ...
    'G9b: the held notch must cover the covariance tracker''s reacquisition lag.');
end


function test_g7_no_regression_on_spatial_tracker(testCase)
% G7: the tone must not disturb the P2 covariance tracker. Runs the same closed
% loop with the waveform layer off and on and compares steady-state oracle gap.
% This is G2's consequence measured where it matters.
td  = testCase.TestData;
scn_off = make_scn(td);
scn_on  = make_scn(td, 'freq', 'constant', 'f_j_offset_hz', 3.1e6);
gap_off = tracker_gap(td, td.sc,        scn_off);
gap_on  = tracker_gap(td, with_fs(td),  scn_on);
fprintf('  G7 steady-state oracle gap: white %.3f dB -> tone %.3f dB (delta %+.3f)\n', ...
    gap_off, gap_on, gap_on - gap_off);
verifyLessThan(testCase, gap_on - gap_off, 0.3, ...
    'G7: the CW jammer must not degrade the spatial tracker.');
end


function test_g10_presence_holds_through_off(testCase)
% G10: with the jammer silent the tracker must HOLD, not chase noise peaks.
%
% Regression for a real bug found on the real array (2026-08-03): the shipped
% presence_snr_db of 10 dB sat INSIDE the noise distribution. A white
% periodogram's own max/median over nfft bins is ~9.5 dB with nothing
% transmitting at all, so 10 dB false-alarmed on 19% of OFF steps and the
% tracker re-locked onto noise, throwing the notch across the band. Measured
% separation on the real array is [12.6 dB OFF max, 18.7 dB ON min]; the tuned
% 15.0 dB sits in that gap. This gate pins the property, not the number.
td  = testCase.TestData;
sc  = with_fs(td);
scn_cfg = struct('id', 'G10', 'motion', 'static', 'power', 'onoff', ...
                 'duty_cycle', 0.5, 'toggle_period_s', 4.0, ...
                 'theta_j_deg', 40.0, 'phi_j_deg', 0.0, ...
                 'freq', 'constant', 'f_j_offset_hz', 3.1e6);
scn = sim_scenario(scn_cfg, td.aj, sc);
cfg = struct('adapt', struct('forgetting_lambda', 0.90, 'diagonal_loading_db', 10, ...
                  'freq', struct('nfft_factor', 4, 'presence_snr_db', 15.0, ...
                                 'smoothing_lambda', 0.7)));
log = closed_loop_run('lcmv', td.stack, [], td.theta_deg, td.phi_deg, scn, td.aj, sc, cfg, []);

on = scn.jammer_on;
fa = mean(log.freq_present(~on));      % false alarms while silent
md = mean(~log.freq_present(on));      % missed detections while transmitting
% Once locked, the held estimate must stay put through every OFF gap.
held = log.f_hat_hz(~on);
held = held(~isnan(held));
spread_hz = max(held) - min(held);
fprintf('  G10 false alarm %.1f%% of OFF steps | missed %.1f%% of ON steps | held estimate spread %.0f Hz\n', ...
    100 * fa, 100 * md, spread_hz);
verifyLessThan(testCase, fa, 0.02, ...
    'G10: the presence test must not fire on a jammer-free spectrum.');
verifyLessThan(testCase, md, 0.02, ...
    'G10: the presence test must still detect a transmitting jammer.');
verifyLessThan(testCase, spread_hz, 5e3, ...
    'G10: the estimate must be HELD through OFF gaps, not re-estimated from noise.');
end


function test_g8_contracts_and_validation(testCase)
% G8: the Mode C contract and every no-silent-defaults guard.
td  = testCase.TestData;
sc  = with_fs(td);
scn = make_scn(td, 'freq', 'constant', 'f_j_offset_hz', 3.1e6);
fr  = adapt_freq_init(struct('freq', struct('nfft_factor', 4, ...
        'presence_snr_db', 15.0, 'smoothing_lambda', 0.7)));

verifyError(testCase, @() adapt_freq_update(fr, struct('sinr_db', 0, 'snapshots', [], 'fs_hz', [])), ...
    'adapt_freq_update:NoSnapshots', 'G8: must reject a Mode S observation.');
verifyError(testCase, @() adapt_freq_update(fr, struct('sinr_db', 0, ...
        'snapshots', ones(8, 128), 'fs_hz', [])), ...
    'adapt_freq_update:NoSampleRate', 'G8: must reject snapshots with no frequency axis.');

% fs_hz set but the scenario carries no carrier track.
verifyError(testCase, @() sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, ...
        make_scn(td), td.aj, sc, 'C'), ...
    'sim_engine_init:MissingCarrier', 'G8: fs_hz without a scenario carrier must error.');

% Carrier outside the captured band would alias.
verifyError(testCase, @() sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, ...
        make_scn(td, 'freq', 'constant', 'f_j_offset_hz', 20e6), td.aj, sc, 'C'), ...
    'sim_engine_init:CarrierOutOfBand', 'G8: an out-of-band carrier must error.');

% Active notch without a frequency axis.
verifyError(testCase, @() sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, ...
        make_scn(td), td.aj, td.sc, 'C', td.notch_cfg), ...
    'sim_engine_init:NotchWithoutWaveform', 'G8: an active notch needs sim.fs_hz.');

% ...but an INERT notch section is fine without one (that is how config.yaml ships).
n_off = td.notch_cfg; n_off.mode = 'off';
st = sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, make_scn(td), td.aj, td.sc, 'C', n_off);
verifyFalse(testCase, st.use_waveform, 'G8: an inert notch must not enable the waveform layer.');

bad = td.notch_cfg; bad.mode = 'notch_it_hard';
verifyError(testCase, @() sim_engine_init(td.stack, [], td.theta_deg, td.phi_deg, ...
        scn, td.aj, sc, 'C', bad), ...
    'sim_engine_init:BadNotchMode', 'G8: an unknown notch mode must error.');

verifyError(testCase, @() adapt_freq_init(struct()), ...
    'adapt_freq_init:MissingKey', 'G8: a missing adapt.freq section must error by name.');
end
