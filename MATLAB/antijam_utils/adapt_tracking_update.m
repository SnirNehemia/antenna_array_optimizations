function [w, state] = adapt_tracking_update(state, obs)
% ADAPT_TRACKING_UPDATE  Mode C update: fold in snapshots, recompute LCMV weights.
%
%   [w, state] = ADAPT_TRACKING_UPDATE(state, obs)
%
%   Consumes obs.snapshots ONLY (Mode C contract; obs.sinr_db is ignored so the
%   tracker never depends on the scalar channel). The K snapshot columns of one
%   step are batch-averaged into a single within-step sample covariance before
%   ONE exponential-forgetting update is applied:
%       R_batch <- mean_i(x_i * x_i'),   i = 1..K
%       R_hat   <- lambda * R_hat + (1 - lambda) * R_batch
%   then w = adapt_lcmv(R_hat, e_s, loading). Note the snapshots contain the
%   desired signal (engine snapshot model), so this is an MPDR-style tracker —
%   diagonal loading is the guard against finite-sample signal self-nulling.
%
%   [P8 fix, 2026-08-01] Previously this applied the (1 - lambda) recursion
%   once PER SNAPSHOT COLUMN (K sequential updates per step) rather than
%   batch-averaging first. Since lambda is calibrated as the inter-step
%   forgetting rate, K sequential within-step recursions silently applied K
%   steps' worth of forgetting inside a single 1-step tick — so raising
%   snapshots_per_step made the estimate WORSE, not better (verified: K
%   16->256 at fixed lambda increased, not decreased, steady-state noise).
%   Batch-averaging first restores the intended meaning of K (more snapshots
%   per step = a better per-step estimate) without changing lambda, so
%   jammer-reacquisition speed after a toggle/jump is unaffected. Measured on
%   the mode_c_demo static-jammer regression: steady-state oracle gap
%   10.4 dB -> 3.4 dB, consecutive-step weight cosine similarity 0.82 -> 0.997
%   at unchanged (loading=10 dB, lambda=0.98, K=16).
%
%   After the closed-form target weights are recomputed, state.mu (see
%   adapt_tracking_init) optionally rate-limits the APPLIED weights toward
%   that target instead of snapping to it directly — see smooth_weights below.
%
%   [P9, 2026-08-02] When state.adaptive_loading is true (see
%   adapt_tracking_init), the loading used above is recomputed each step as
%   loading_factor * sqrt(sig_power_hat * noise_floor_hat) instead of staying
%   fixed — see adapt_tracking_init.m's header note for why the geometric mean
%   of the measured signal power and noise floor, not the noise floor alone,
%   is the right reference (noise_floor_estimate below computes the latter).
%
%   [P11, 2026-08-31] sig_power_hat is now the CAPON estimate
%   (capon_power_estimate below), not the Bartlett quotient P9 shipped with —
%   that one read the jammer rather than the signal and over-loaded the
%   beamformer into a null-starved collapse. The result is additionally
%   floored at the fixed diagonal_loading_db (state.loading_floor). See
%   capon_power_estimate and adapt_tracking_init's header for the evidence.
%
%   Inputs:
%       state : struct from adapt_tracking_init.
%       obs   : observation struct from sim_engine_step (Mode C).
%
%   Outputs:
%       w     : (N_el x 1) complex weights to apply next step.
%       state : updated state (R_hat, w).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P2, P9].

if isempty(obs.snapshots)
    error('adapt_tracking_update:NoSnapshots', ...
        'obs.snapshots is empty — the covariance tracker requires Mode C.');
end

X = obs.snapshots;
R_batch = (X * X') / size(X, 2);
state.R_hat = state.lambda * state.R_hat + (1 - state.lambda) * R_batch;

if state.adaptive_loading
    sigma_n_hat_raw = noise_floor_estimate(state.R_hat, state.n_sig);
    state.noise_floor_hat = state.lambda * state.noise_floor_hat + ...
        (1 - state.lambda) * sigma_n_hat_raw;
    % [P11, 2026-08-31] Capon (MVDR) estimate of the desired-signal power, NOT
    % the Bartlett/conventional one this used to be. See the estimator note in
    % adapt_tracking_init's header for the measurement that forced the change.
    sig_power_raw = capon_power_estimate(state.R_hat, state.e_s, ...
        state.sig_power_hat);
    state.sig_power_hat = state.lambda * state.sig_power_hat + ...
        (1 - state.lambda) * sig_power_raw;
    % [P11] Floored at the fixed diagonal_loading_db (state.loading_floor) —
    % see adapt_tracking_init's header: unfloored, this collapses at the
    % signal-limited edge where sig_power_hat stops separating from the noise.
    state.loading = max(state.loading_floor, state.loading_factor * ...
        sqrt(state.sig_power_hat * state.noise_floor_hat));
end

w_target = adapt_lcmv(state.R_hat, state.e_s, state.loading);
w = smooth_weights(state.w, w_target, state.mu);
state.w = w;
end


function sigma_n_sq = noise_floor_estimate(R_hat, n_sig)
% NOISE_FLOOR_ESTIMATE  Median of R_hat's noise-subspace eigenvalues.
%   Eigendecomposes the (forced-Hermitian) covariance, sorts descending, and
%   takes the median (not the mean — robust to a thin noise subspace, e.g. a
%   6-element array leaves only N_el - n_sig ~ 2 eigenvalues to summarize) of
%   the bottom N_el - n_sig eigenvalues. Same split as adapt_music_doa.m's
%   signal/noise subspace (duplicated, not shared — see adapt_tracking_init.m
%   header note).
n_el = size(R_hat, 1);
lam  = eig((R_hat + R_hat') / 2);
lam  = sort(real(lam), 'descend');
sigma_n_sq = median(lam((n_sig + 1):n_el));
end


% ────────────────────────── HELPERS ───────────────────────────────

function w = smooth_weights(w_prev, w_target, mu)
% SMOOTH_WEIGHTS  Rate-limit the applied weights toward a closed-form target.
%   mu = 1 (default, feature off) returns w_target unchanged. mu < 1 blends
%   w_prev and w_target after aligning their global phase — MVDR/max-SINR
%   solutions are only defined up to an arbitrary unit-modulus phase, so a
%   naive blend could destructively interfere between two representations of
%   the "same" beam and corrupt the result.
if mu >= 1
    w = w_target;
    return
end
ph = w_prev' * w_target;
ph = ph / max(abs(ph), eps);
w_target_aligned = w_target * conj(ph) / abs(ph);
w = (1 - mu) * w_prev + mu * w_target_aligned;
w = w / norm(w);
end


function p = capon_power_estimate(R_hat, e_s, p_prev)
% CAPON_POWER_ESTIMATE  Interference-nulling estimate of the desired-signal power.
%
%   p = trace(inv(e_s' * inv(R_hat) * e_s)) / n_comp
%
%   This is the Capon (MVDR) spectral estimator evaluated at the steering
%   direction: it is the output power of the minimum-variance distortionless
%   beamformer pointed at e_s, so any source NOT at e_s is nulled before the
%   power is read.
%
%   [P11, 2026-08-31] It replaces the Bartlett (conventional-beamformer)
%   estimator trace(e_s' * R_hat * e_s) / trace(e_s' * e_s) that P9 originally
%   used. Bartlett applies the QUIESCENT beam and reads whatever that beam
%   collects, which includes the jammer through the sidelobes — so it does not
%   estimate the desired signal at all once the jammer is strong. Measured on
%   the real array (patchs_with_monopoles, total-pol, sigma_s_db = 0, i.e. a
%   true desired power of 1) the Bartlett estimate read +42.5 / +52.0 / +61.9
%   dB as jn_ratio_db went 10 / 20 / 30 — it tracked the JAMMER, one-for-one.
%   The loading computed from it therefore climbed to ~17.7 dB at J/N = 20,
%   against a jammer eigenvalue of only ~27.8 dB (J/N + 10*log10(N_el)), which
%   is exactly the over-loading that config.yaml's P2 note warns starves the
%   null — and the P11 campaign measured the consequence as a total collapse
%   (0% availability) for sigma_s_db <= 4 at high J/S. The Capon estimate at
%   the same three points reads 1.09 / 1.13 / 1.14: invariant to the jammer,
%   which is the whole requirement.
%
%   p_prev is used only if R_hat is too ill-conditioned to invert reliably
%   (it should not be — R_hat is seeded at eye(.) and stays positive definite
%   under the convex forgetting update — but an overnight campaign should
%   degrade rather than abort).
n_comp = size(e_s, 2);
R_sym  = (R_hat + R_hat') / 2;
if ~isfinite(rcond(R_sym)) || rcond(R_sym) < 1e-12
    p = p_prev;
    return
end
M = e_s' * (R_sym \ e_s);
M = (M + M') / 2;
if ~isfinite(rcond(M)) || rcond(M) < 1e-12
    p = p_prev;
    return
end
p = real(trace(M \ eye(n_comp))) / n_comp;
if ~isfinite(p) || p <= 0
    p = p_prev;
end
end
