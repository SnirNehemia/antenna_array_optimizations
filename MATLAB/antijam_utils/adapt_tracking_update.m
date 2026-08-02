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
%   Inputs:
%       state : struct from adapt_tracking_init.
%       obs   : observation struct from sim_engine_step (Mode C).
%
%   Outputs:
%       w     : (N_el x 1) complex weights to apply next step.
%       state : updated state (R_hat, w).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P2].

if isempty(obs.snapshots)
    error('adapt_tracking_update:NoSnapshots', ...
        'obs.snapshots is empty — the covariance tracker requires Mode C.');
end

X = obs.snapshots;
R_batch = (X * X') / size(X, 2);
state.R_hat = state.lambda * state.R_hat + (1 - state.lambda) * R_batch;

w_target = adapt_lcmv(state.R_hat, state.e_s, state.loading);
w = smooth_weights(state.w, w_target, state.mu);
state.w = w;
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
