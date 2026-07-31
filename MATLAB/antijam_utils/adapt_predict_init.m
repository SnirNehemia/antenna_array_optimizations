function state = adapt_predict_init(adapt_config, aj_config, e_s, E1, E2, theta_deg, phi_deg, n_elements)
% ADAPT_PREDICT_INIT  State for the Mode C predictive (anticipatory) nuller.
%
%   state = ADAPT_PREDICT_INIT(adapt_config, aj_config, e_s, E1, E2, ...
%                              theta_deg, phi_deg, n_elements)
%
%   Initializes the predictive Mode C algorithm: a recursive covariance
%   estimate (same forgetting buffer as adapt_tracking), a ring buffer of the
%   MUSIC jammer-presence indicator, and periodogram state used to learn the
%   jammer's on/off duty period and pre-form the null before it turns back on.
%   R_hat starts at the identity (quiescent MVDR first weights), matching
%   adapt_tracking_init. Retains the grid (E1/E2/theta/phi) so each update can
%   run adapt_music_doa and build the null steering column at the predicted
%   angle. Consumes obs.snapshots ONLY (Mode C contract).
%
%   Inputs:
%       adapt_config : struct. Required: forgetting_lambda, diagonal_loading_db,
%                      and the nested predict block with fields
%                      presence_gap_db, buffer_len, min_periods, lead_steps,
%                      doa_stride.
%       aj_config    : antijam config section — supplies theta_s_deg, phi_s_deg,
%                      guard_deg for the MUSIC guard mask (doa cfg).
%       e_s          : (N_el x n_c) steering column(s) toward theta_s.
%       E1           : (N_el x N_theta x N_phi) primary far-field stack.
%       E2           : (N_el x N_theta x N_phi) secondary component, or [].
%       theta_deg    : (1 x N_theta) elevation grid [deg].
%       phi_deg      : (1 x N_phi) azimuth grid [deg].
%       n_elements   : N_el.
%
%   Outputs:
%       state : struct with fields
%           R_hat      : (N_el x N_el) running covariance estimate.
%           lambda     : forgetting factor (per snapshot column).
%           loading    : linear diagonal loading (from diagonal_loading_db).
%           e_s        : constraint column(s).
%           E1,E2,theta_deg,phi_deg : grid refs (MUSIC + null steering).
%           doa_cfg    : struct passed to adapt_music_doa each step.
%           presence   : (1 x k) growing 0/1 record of the presence indicator
%                        (the periodogram analyses its last buffer_len samples).
%           buffer_len : analysis-window length for the periodogram [steps].
%           min_periods: cycles to observe before trusting the period estimate.
%           lead_steps : pre-form the null this many steps before predicted ON.
%           period_est : detected duty period [steps] (NaN until learned).
%           last_doa   : last-known jammer (theta,phi) while present [deg].
%           k          : step counter.
%           w          : (N_el x 1) current weights (initial quiescent MVDR).
%           diag       : struct of per-step diagnostics appended for the report
%                        (theta_j_deg, phi_j_deg, present, predicted_on,
%                        period_est, pspec) — see plot_doa_waterfall.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

REQUIRED = {'forgetting_lambda', 'diagonal_loading_db'};
for i = 1:numel(REQUIRED)
    if ~isfield(adapt_config, REQUIRED{i}) || isempty(adapt_config.(REQUIRED{i}))
        error('adapt_predict_init:MissingKey', ...
            'Missing required adapt config key: ''%s''.', REQUIRED{i});
    end
end
if ~isfield(adapt_config, 'predict') || isempty(adapt_config.predict)
    error('adapt_predict_init:MissingKey', ...
        'Missing required adapt config block: ''predict''.');
end
p = adapt_config.predict;
PREQ = {'presence_gap_db', 'buffer_len', 'min_periods', 'lead_steps', 'doa_stride'};
for i = 1:numel(PREQ)
    if ~isfield(p, PREQ{i}) || isempty(p.(PREQ{i}))
        error('adapt_predict_init:MissingKey', ...
            'Missing required adapt.predict key: ''%s''.', PREQ{i});
    end
end

state = struct();
state.R_hat   = eye(n_elements);                  % quiescent start (unit noise)
state.lambda  = adapt_config.forgetting_lambda;
state.loading = 10^(adapt_config.diagonal_loading_db / 10);
state.e_s     = e_s;

% Grid references — needed each step to run MUSIC and to build the steering
% column at the (predicted) jammer angle for the hard null constraint.
state.E1        = E1;
state.E2        = E2;
state.theta_deg = theta_deg;
state.phi_deg   = phi_deg;
state.n_el      = n_elements;

% Config passed straight to adapt_music_doa each step.
state.doa_cfg = struct('theta_s_deg', aj_config.theta_s_deg, ...
    'phi_s_deg', aj_config.phi_s_deg, 'guard_deg', aj_config.guard_deg, ...
    'presence_gap_db', p.presence_gap_db, 'doa_stride', p.doa_stride, ...
    'return_pspec', false);

% Presence history + on/off period-detection state. The periodogram analyses
% the most recent buffer_len samples of this growing 0/1 record.
state.buffer_len  = p.buffer_len;
state.presence    = [];                            % 1 = jammer detected, per step
state.min_periods = p.min_periods;
state.lead_steps  = p.lead_steps;
state.period_est  = NaN;                           % [steps], NaN until learned
state.duty_est    = NaN;                           % ON fraction, NaN until learned
state.last_on_k   = NaN;                           % step index of last 0->1 turn-on
state.last_doa    = struct('theta_deg', NaN, 'phi_deg', NaN, 'idx', NaN);
state.k           = 0;

% Per-step diagnostics for closed_loop_run / plot_doa_waterfall.
state.last = struct('theta_j_deg', NaN, 'phi_j_deg', NaN, 'present', false, ...
    'predicted_on', false, 'period_est', NaN);

% Initial weights: quiescent MVDR (no jammer seen yet).
state.w = adapt_lcmv_null(state.R_hat, e_s, [], state.loading);
end
