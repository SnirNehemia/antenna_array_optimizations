function state = adapt_tracking_init(adapt_config, e_s, n_elements)
% ADAPT_TRACKING_INIT  State for the Mode C recursive-covariance LCMV tracker.
%
%   state = ADAPT_TRACKING_INIT(adapt_config, e_s, n_elements)
%
%   Initializes R_hat = sigma_n^2 * I (= identity under the engine's unit noise
%   floor) so the first weights are the quiescent MVDR solution.
%
%   Inputs:
%       adapt_config : struct. Required: forgetting_lambda, diagonal_loading_db.
%       e_s          : (N_el x n_c) steering column(s) toward theta_s.
%       n_elements   : N_el.
%
%   Outputs:
%       state : struct with fields
%           R_hat   : (N_el x N_el) running covariance estimate.
%           lambda  : forgetting factor (applied per snapshot column).
%           loading : linear diagonal loading (from diagonal_loading_db,
%                     relative to the sigma_n^2 = 1 noise floor).
%           mu      : weight-vector smoothing factor, see
%                     adapt_config.weight_smoothing_mu below.
%           e_s     : constraint column(s).
%           w       : (N_el x 1) current weights (initial quiescent solution).
%
%   [P8, 2026-08-01] Optional weight_smoothing_mu (OPT-IN, not in REQUIRED):
%   adapt_lcmv recomputes a full closed-form solution from R_hat every step —
%   a "snap to the current estimate," not an iterative climb — so even a
%   well-converged R_hat can still move the applied w noticeably between
%   consecutive steps when the SINR optimum is sharp (see docs/notes.md [P8]
%   2026-08-01 entry). mu < 1 rate-limits the APPLIED weights toward the
%   target solution (adapt_tracking_update: w <- (1-mu)*w + mu*w_target,
%   phase-aligned first), trading some mean oracle gap and reacquisition speed
%   for a visibly gradual approach. Absent/empty -> mu = 1 (no smoothing,
%   original snap-to-target behavior) — this is a deliberate default (the
%   feature is opt-in), not the "no silent defaults" case, which applies to
%   forgetting_lambda/diagonal_loading_db above.
%
%   [P9, 2026-08-02] Optional loading_factor_db (OPT-IN, not in REQUIRED):
%   diagonal_loading_db fixes the loading as an ASSUMED sigma_n^2 = 1 dB
%   offset; that assumption breaks whenever sigma_s_db/jn_ratio_db land far
%   from the regime it was swept against (see docs/notes.md [P9] entry — the
%   patchs_with_monopoles/sigma_s_db=30 regression, where the desired signal
%   sits 20 dB ABOVE the jammer instead of well below it). A first attempt at
%   scaling loading off the estimated NOISE floor alone did not fix the
%   regression (verified empirically: noise floor is fixed at sigma_n^2=1 by
%   construction regardless of regime, so that gave the same numeric loading
%   as the old fixed value). The corrected quantity is the GEOMETRIC MEAN of
%   the noise floor and the desired-signal power actually measured in R_hat:
%       loading = loading_factor * sqrt(sig_power_hat * noise_floor_hat)
%   sig_power_hat is the Rayleigh quotient of R_hat at the KNOWN e_s direction
%   (trace(e_s' R_hat e_s) / trace(e_s' e_s)) — cheap, and sidesteps any
%   signal/jammer subspace-rank assumption since e_s (unlike the jammer angle)
%   is not unknown. noise_floor_hat is the median of the bottom
%   N_el - 2*n_comp "noise" eigenvalues of R_hat (n_comp = size(e_s,2); the
%   top 2*n_comp span the desired signal + jammer, matching adapt_music_doa.m's
%   model-order-2 convention — duplicated here, not shared, so this module has
%   no dependency on the P8 MUSIC file). Both estimates are smoothed with the
%   same forgetting_lambda EMA as R_hat itself (state.sig_power_hat,
%   state.noise_floor_hat) so per-step eigenvalue/Rayleigh-quotient jitter
%   doesn't reintroduce jitter into the loading. Empirically verified
%   (mode_c_demo, 2026-08-02): reproduces the P2-tuned ~10 linear in the
%   original weak-signal regime (sigma_s_db=3) with loading_factor_db=0, and
%   lands at ~205 in the broken sigma_s_db=30 regime — within the empirical
%   sweet spot found by a brute-force fixed-loading sweep on that scenario
%   (best ~100-300, oracle gap ~8.2 dB vs 13.6 dB at the old fixed 10). Absent
%   loading_factor_db -> adaptive loading OFF (state.adaptive_loading = false),
%   the original fixed-diagonal_loading_db behavior — opt-in, like
%   weight_smoothing_mu above; diagonal_loading_db stays REQUIRED regardless
%   since it is still the value used whenever loading_factor_db is absent.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P2, P9].

REQUIRED = {'forgetting_lambda', 'diagonal_loading_db'};
for i = 1:numel(REQUIRED)
    if ~isfield(adapt_config, REQUIRED{i}) || isempty(adapt_config.(REQUIRED{i}))
        error('adapt_tracking_init:MissingKey', ...
            'Missing required adapt config key: ''%s''.', REQUIRED{i});
    end
end

state = struct();
state.R_hat   = eye(n_elements);
state.lambda  = adapt_config.forgetting_lambda;
state.loading = 10^(adapt_config.diagonal_loading_db / 10);
if isfield(adapt_config, 'weight_smoothing_mu') && ~isempty(adapt_config.weight_smoothing_mu)
    state.mu = adapt_config.weight_smoothing_mu;
else
    state.mu = 1.0;                                    % opt-in feature: off by default
end
state.e_s     = e_s;

% [P9] Data-driven loading (opt-in — see header note above).
if isfield(adapt_config, 'loading_factor_db') && ~isempty(adapt_config.loading_factor_db)
    n_comp = size(e_s, 2);
    n_sig  = 2 * n_comp;                   % desired signal + 1 jammer (locked scope)
    if n_sig >= n_elements
        error('adapt_tracking_init:TooFewElements', ...
            'Adaptive loading needs N_el > 2*n_comp (= %d); got N_el = %d.', ...
            n_sig, n_elements);
    end
    state.adaptive_loading = true;
    state.n_sig            = n_sig;
    state.loading_factor   = 10^(adapt_config.loading_factor_db / 10);
    state.noise_floor_hat  = 1.0;           % matches R_hat = eye(.) at k=0
    state.sig_power_hat    = 1.0;           % Rayleigh quotient of eye(.) at e_s
    state.loading          = state.loading_factor * ...
        sqrt(state.sig_power_hat * state.noise_floor_hat);
else
    state.adaptive_loading = false;
end

state.w       = adapt_lcmv(state.R_hat, e_s, state.loading);
end
