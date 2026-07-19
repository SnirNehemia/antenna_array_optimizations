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
%           e_s     : constraint column(s).
%           w       : (N_el x 1) current weights (initial quiescent solution).
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
state.e_s     = e_s;
state.w       = adapt_lcmv(state.R_hat, e_s, state.loading);
end
