function [w, state] = adapt_predict_update(state, obs)
% ADAPT_PREDICT_UPDATE  Mode C update: MUSIC DoA + on/off prediction -> pre-null.
%
%   [w, state] = ADAPT_PREDICT_UPDATE(state, obs)
%
%   Consumes obs.snapshots ONLY (Mode C contract; obs.sinr_db is ignored so the
%   algorithm never depends on the scalar channel). Per step:
%     1. Batch-average the K snapshot columns into a within-step sample
%        covariance, then fold it into R_hat with ONE exponential-forgetting
%        update (identical accumulation to adapt_tracking_update; see the
%        [P8 fix, 2026-08-01] note there for why per-column recursion was
%        wrong).
%     2. doa = adapt_music_doa(R_hat, ...): jammer (theta,phi), presence, power.
%     3. Push doa.present into the presence ring buffer; if present, refresh
%        state.last_doa.
%     4. Periodogram (base-MATLAB fft) of the mean-removed presence buffer once
%        >= min_periods cycles are available: dominant non-DC bin -> period_est;
%        bin phase -> predicted next turn-on. predicted_on is true when the
%        jammer is expected ON within lead_steps.
%     5. Choose the null constraint and recompute weights:
%          - jammer present  -> null at the MEASURED doa (theta_j,phi_j), R = R_hat;
%          - predicted ON    -> pre-null at last_doa, R = eye (quiescent), so the
%                               null is formed before any energy returns;
%          - otherwise       -> no null (plain MVDR / distortionless beam).
%        via w = adapt_lcmv_null(R, e_s, e_null, loading).
%     6. Append this step's diagnostics to state.diag for plot_doa_waterfall.
%
%   state.mu (see adapt_predict_init) optionally rate-limits the APPLIED
%   weights toward whichever w_target the branch above selected, instead of
%   snapping to it directly — see smooth_weights below (shared logic with
%   adapt_tracking_update.m).
%
%   S6 note: a single ON burst has no periodic line — the periodogram stays
%   unlearned and the null is held from last_doa (persist-last fallback). S5's
%   simultaneous drift makes the held-static angle stale over the OFF gap; fully
%   closing S5 needs the CV-Kalman drift predictor (P8 follow-up).
%
%   Inputs:
%       state : struct from adapt_predict_init (advanced state returned).
%       obs   : observation struct from sim_engine_step (Mode C; snapshots set).
%
%   Outputs:
%       w     : (N_el x 1) complex weights to apply next step.
%       state : updated state (R_hat, ring buffer, period_est, last_doa, k, w,
%               diag).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

if isempty(obs.snapshots)
    error('adapt_predict_update:NoSnapshots', ...
        'obs.snapshots is empty — the predictive tracker requires Mode C.');
end
state.k = state.k + 1;
k = state.k;

% ── 1. Fold snapshots into R_hat (same forgetting as adapt_tracking) ──
X = obs.snapshots;
R_batch = (X * X') / size(X, 2);
state.R_hat = state.lambda * state.R_hat + (1 - state.lambda) * R_batch;

% [P9] Data-driven loading (opt-in — see adapt_tracking_init.m/
% adapt_tracking_update.m for the full geometric-mean rationale; duplicated
% here, not shared).
if state.adaptive_loading
    sigma_n_hat_raw = noise_floor_estimate(state.R_hat, state.n_sig);
    state.noise_floor_hat = state.lambda * state.noise_floor_hat + ...
        (1 - state.lambda) * sigma_n_hat_raw;
    sig_power_raw = real(trace(state.e_s' * state.R_hat * state.e_s)) / ...
        real(trace(state.e_s' * state.e_s));
    state.sig_power_hat = state.lambda * state.sig_power_hat + ...
        (1 - state.lambda) * sig_power_raw;
    state.loading = state.loading_factor * ...
        sqrt(state.sig_power_hat * state.noise_floor_hat);
end

% ── 2. MUSIC: explicit jammer angle + presence from this R_hat ─────
doa = adapt_music_doa(state.R_hat, state.E1, state.E2, ...
    state.theta_deg, state.phi_deg, state.doa_cfg);

% ── 3. Update presence history + last-known jammer direction ───────
prev_present = ~isempty(state.presence) && state.presence(end) > 0;
state.presence(end + 1) = double(doa.present);
if doa.present && ~prev_present
    state.last_on_k = k;                          % 0 -> 1 turn-on this step
end
if doa.present
    [it, ip] = nearest_index_2d(state.theta_deg, state.phi_deg, ...
        doa.theta_j_deg, doa.phi_j_deg);
    state.last_doa = struct('theta_deg', doa.theta_j_deg, ...
        'phi_deg', doa.phi_j_deg, ...
        'idx', (ip - 1) * numel(state.theta_deg) + it);
end

% ── 4. Periodogram: learn the on/off toggle period + duty cycle ────
% The customer's "spectrogram of repeating patterns": FFT the recent presence
% record, take the dominant non-DC line as the toggle frequency.
win = state.presence(max(1, end - state.buffer_len + 1):end);
[period_est, duty_est, trusted] = detect_period(win, state.min_periods);
if trusted
    state.period_est = period_est;
    state.duty_est   = duty_est;
end

% ── 5. Predict whether the jammer will be ON lead_steps from now ───
% Present now => obviously on. Otherwise, if the period is trusted, fire the
% pre-null in the lead_steps window BEFORE the next projected turn-on (the next
% turn-on is last_on_k + one period; steps_to_next counts down to it).
predicted_on = doa.present;
if ~doa.present && trusted && ~isnan(state.last_on_k)
    ph = mod(k - state.last_on_k, state.period_est);   % steps since last on-start
    steps_to_next = state.period_est - ph;             % steps to next on-start
    predicted_on = steps_to_next <= state.lead_steps;
end

% ── 6. Choose the null constraint and recompute the weights ────────
% Design invariant: while the jammer is visible, let the measured R_hat form
% the null (exactly the reactive LCMV tracker) — this makes 'predict' never
% worse than 'lcmv' during ON, and avoids mis-steering on a transient MUSIC
% estimate. The explicit hard null is used ONLY to pre-form the null during a
% predicted-imminent OFF window, where R_hat carries no jammer energy.
if doa.present
    % Jammer visible: reactive MPDR null from the measured covariance.
    w_target = adapt_lcmv_null(state.R_hat, state.e_s, [], state.loading);
elseif predicted_on && ~isnan(state.last_doa.idx)
    % About to return: pre-null the LAST-KNOWN direction from a quiescent R
    % (identity) — the hard constraint forms the null before any energy is back.
    e_null = steer_col(state, state.last_doa.idx);
    w_target = adapt_lcmv_null(eye(state.n_el), state.e_s, e_null, state.loading);
else
    % Jammer absent and no imminent return: restore the quiescent full-gain beam.
    w_target = adapt_lcmv_null(eye(state.n_el), state.e_s, [], state.loading);
end
w = smooth_weights(state.w, w_target, state.mu);
state.w = w;

% ── 7. Record per-step diagnostics ─────────────────────────────────
state.last = struct('theta_j_deg', doa.theta_j_deg, 'phi_j_deg', doa.phi_j_deg, ...
    'present', doa.present, 'predicted_on', predicted_on, ...
    'period_est', state.period_est);
end


% ────────────────────────── HELPERS ───────────────────────────────

function e = steer_col(state, idx)
% Steering column(s) at flat grid index idx (both polarization components).
e = state.E1(:, idx);
if ~isempty(state.E2)
    e = [e, state.E2(:, idx)];
end
end


function [period, duty, trusted] = detect_period(p, min_periods)
% On/off period + duty from the presence record via a base-MATLAB periodogram.
%   period  : dominant cycle length [steps] (NaN if not learned).
%   duty    : fraction of samples ON (the projected ON-window width).
%   trusted : true once a clear line is found AND >= min_periods cycles seen.
period = NaN; duty = mean(p); trusted = false;
n = numel(p);
if n < 8 || all(p == p(1))                        % too short / no toggling yet
    return
end
P = abs(fft(p - mean(p))).^2;                      % periodogram (power spectrum)
half = floor(n / 2);
[pk, bin] = max(P(2:half + 1));                    % dominant non-DC bin
bin = bin + 1;                                     % undo the 2: offset
if pk <= 3 * mean(P(2:half + 1))                   % peak not prominent -> no line
    return
end
% Parabolic (sub-bin) interpolation of the peak for an accurate period: the FFT
% bin grid is coarse (period = n/bin), so the raw bin can be off by several
% steps — enough to mistime the pre-null.
delta = 0;
if bin > 2 && bin < half + 1
    a = P(bin - 1); b = P(bin); c = P(bin + 1);
    delta = 0.5 * (a - c) / (a - 2 * b + c + eps); % in [-0.5, 0.5]
end
period = n / (bin - 1 + delta);                    % bin frequency -> period [steps]
trusted = n >= min_periods * period;               % enough cycles observed
end


function sigma_n_sq = noise_floor_estimate(R_hat, n_sig)
% NOISE_FLOOR_ESTIMATE  Median of R_hat's noise-subspace eigenvalues.
%   See adapt_tracking_update.m's copy of this helper for the full rationale
%   (duplicated, not shared, per the P9 design decision).
n_el = size(R_hat, 1);
lam  = eig((R_hat + R_hat') / 2);
lam  = sort(real(lam), 'descend');
sigma_n_sq = median(lam((n_sig + 1):n_el));
end


function w = smooth_weights(w_prev, w_target, mu)
% SMOOTH_WEIGHTS  Rate-limit the applied weights toward a closed-form target.
%   mu = 1 (default, feature off) returns w_target unchanged. mu < 1 blends
%   w_prev and w_target after aligning their global phase — MVDR/max-SINR/
%   null-constrained solutions are only defined up to an arbitrary
%   unit-modulus phase, so a naive blend could destructively interfere
%   between two representations of the "same" beam and corrupt the result.
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
