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
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P2].

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
state.w       = adapt_lcmv(state.R_hat, e_s, state.loading);
end
