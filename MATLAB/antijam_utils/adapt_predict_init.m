function state = adapt_predict_init(adapt_config, aj_config, e_s, E1, E2, theta_deg, phi_deg, n_elements, sim_config)
% ADAPT_PREDICT_INIT  State for the Mode C predictive (anticipatory) nuller.
%
%   state = ADAPT_PREDICT_INIT(adapt_config, aj_config, e_s, E1, E2, ...
%                              theta_deg, phi_deg, n_elements, sim_config)
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
%       sim_config   : sim config section. Read only for dt_s, and only when
%                      the optional adapt.predict.cv block is present ([P12b]).
%
%   Outputs:
%       state : struct with fields
%           R_hat      : (N_el x N_el) running covariance estimate.
%           lambda     : forgetting factor (per snapshot column).
%           loading    : linear diagonal loading (from diagonal_loading_db).
%           mu         : weight-vector smoothing factor. Optional
%                        adapt_config.weight_smoothing_mu; absent -> 1 (off).
%                        See adapt_tracking_init.m for the full rationale —
%                        applies uniformly to whichever of the three w_target
%                        branches adapt_predict_update selects each step.
%           loading_floor : [P11] the fixed diagonal_loading_db value, which
%                        the data-driven loading may never go below.
%           adaptive_loading, n_sig, loading_factor, noise_floor_hat :
%                        [P9] optional data-driven loading, opt-in via
%                        adapt_config.loading_factor_db — see
%                        adapt_tracking_init.m for the full rationale
%                        (duplicated here, not shared).
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
if isfield(adapt_config, 'weight_smoothing_mu') && ~isempty(adapt_config.weight_smoothing_mu)
    state.mu = adapt_config.weight_smoothing_mu;
else
    state.mu = 1.0;                                    % opt-in feature: off by default
end
state.e_s     = e_s;

% [P9] Data-driven loading (opt-in — see adapt_tracking_init.m for rationale).
if isfield(adapt_config, 'loading_factor_db') && ~isempty(adapt_config.loading_factor_db)
    n_comp = size(e_s, 2);
    n_sig  = 2 * n_comp;                   % desired signal + 1 jammer (locked scope)
    if n_sig >= n_elements
        % [O, 2026-09-07] DEGRADE, do not throw. This precondition used to be a
        % hard error, which made the whole Mode C stack unusable on small
        % apertures: `Dipole` (1 el) and `patch_back2back` (2 el, dual-pol) could
        % not run even the REACTIVE tracker, because this opt-in feature aborted
        % initialization. P9 measured adaptive loading at +0.00 pp against fixed
        % on all 90 campaign cells, so refusing to run rather than falling back
        % traded a real capability for no measured benefit. The configured
        % diagonal_loading_db is used instead, and the fallback is announced --
        % it is a warning, not a silent default (CLAUDE.md rule 4).
        warning('adapt_predict_init:AdaptiveLoadingInfeasible', ...
            ['Adaptive loading needs N_el > 2*n_comp (= %d) but this array has ' ...
             'N_el = %d; falling back to the configured diagonal_loading_db ' ...
             '(%.1f dB) for this run.'], n_sig, n_elements, ...
            adapt_config.diagonal_loading_db);
        state.adaptive_loading = false;
    else
        state.adaptive_loading = true;
        state.n_sig            = n_sig;
        state.loading_factor   = 10^(adapt_config.loading_factor_db / 10);
        state.noise_floor_hat  = 1.0;           % matches R_hat = eye(.) at k=0
        state.sig_power_hat    = 1.0;           % Rayleigh quotient of eye(.) at e_s
        % [P11] Floor the data-driven value at the fixed diagonal_loading_db —
        % see adapt_tracking_init's header for the campaign evidence.
        state.loading_floor    = state.loading;
        state.loading          = max(state.loading_floor, state.loading_factor * ...
            sqrt(state.sig_power_hat * state.noise_floor_hat));
    end
else
    state.adaptive_loading = false;
end

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

% [O, 2026-09-07] On/off detection repair. OPT-IN: with the adapt.predict.onoff
% block absent this is inert and the pre-P12c behaviour is byte-identical.
%
% WHY (measured, patchs_with_monopoles, J/N 25 dB, duty 0.5): the anticipatory
% pre-null branch -- the entire reason P8 exists -- fires on AT MOST 1.3% of
% steps, and never at all outside a narrow 10-15 s toggle band. Two independent
% causes, both structural rather than tuning:
%
%   A. SLOW toggling (period >= 20 s): buffer_len = 1024 steps = 51.2 s and
%      min_periods = 3 requires three cycles inside the analysis window, so a
%      period above ~17 s can NEVER be trusted. Arithmetic, not tuning.
%   B. FAST toggling (period <= 5 s): presence saturates at 100%. The
%      covariance horizon is 1/(1-lambda) = 10 steps = 0.5 s but the OFF phase
%      is only 1-2.5 s, so R_hat still carries jammer energy and the eigengap
%      never collapses. The presence signal the periodogram consumes is itself
%      low-pass filtered by lambda. Constant presence -> no spectral line.
%
% The repair is therefore two-part: size the analysis window from the longest
% period worth detecting (A), and derive presence from a SECOND, much faster
% covariance so on/off transitions are visible at all (B). A single forgetting
% factor cannot serve both the beamformer and the presence detector -- the
% beamformer wants a long memory, the detector wants a short one.
if isfield(p, 'onoff') && ~isempty(p.onoff)
    OREQ = {'fast_lambda', 'max_period_s', 'lead_frac'};
    for i = 1:numel(OREQ)
        if ~isfield(p.onoff, OREQ{i}) || isempty(p.onoff.(OREQ{i}))
            error('adapt_predict_init:MissingKey', ...
                'Missing required adapt.predict.onoff key: ''%s''.', OREQ{i});
        end
    end
    if nargin < 9 || isempty(sim_config) || ~isfield(sim_config, 'dt_s')
        error('adapt_predict_init:MissingSimConfig', ...
            'adapt.predict.onoff needs the sim config (9th argument) for dt_s.');
    end
    if p.onoff.fast_lambda >= adapt_config.forgetting_lambda
        error('adapt_predict_init:BadFastLambda', ...
            ['adapt.predict.onoff.fast_lambda (%.3f) must be SHORTER-memory than ' ...
             'forgetting_lambda (%.3f); a detector slower than the beamformer ' ...
             'cannot resolve transitions the beamformer already smooths.'], ...
            p.onoff.fast_lambda, adapt_config.forgetting_lambda);
    end
    state.onoff_enabled = true;
    state.fast_lambda   = p.onoff.fast_lambda;
    state.lead_frac     = p.onoff.lead_frac;
    % Cap the lead at lead_cap_horizons covariance horizons. 1/(1-lambda) is
    % how long R_hat needs to re-acquire the jammer after it returns, so
    % leading further than a small multiple of that pre-nulls into empty space.
    if isfield(p.onoff, 'lead_cap_horizons') && ~isempty(p.onoff.lead_cap_horizons)
        cap_h = p.onoff.lead_cap_horizons;
    else
        cap_h = 2.0;
    end
    state.lead_cap_steps = cap_h / max(1 - adapt_config.forgetting_lambda, eps);
    % How long the OFF window must be, in covariance horizons, before dropping
    % the null is worth the re-acquisition cost at the next turn-on.
    %
    % DEFAULT 0 = always release, which is the un-gated behaviour. The gate is
    % available but is NOT recommended on the evidence: it removes a small tail
    % of regressions at the fastest toggle (a 4 s period cell goes 89.4 -> 90.0
    % instead of 89.4 -> 80.1) but costs far more where the repair pays most --
    % a hard 10 s cell collapses from 62.2 back to 33.0, because holding the
    % null through OFF forfeits the quiescent gain that made the repair
    % worthwhile. Net over the campaign the un-gated form is +7.5 pp; the gate
    % trades that away to tidy the tail. Kept as a knob, defaulted off.
    if isfield(p.onoff, 'release_min_horizons') && ~isempty(p.onoff.release_min_horizons)
        state.release_min_horizons = p.onoff.release_min_horizons;
    else
        state.release_min_horizons = 0.0;
    end
    state.R_fast        = eye(n_elements);
    % (A) size the window so min_periods cycles of the LONGEST period fit.
    state.buffer_len = max(p.buffer_len, ...
        ceil(p.min_periods * p.onoff.max_period_s / sim_config.dt_s));
else
    state.onoff_enabled  = false;
    state.fast_lambda    = NaN;
    state.lead_frac      = NaN;
    state.lead_cap_steps = NaN;
    state.release_min_horizons = NaN;
    state.R_fast         = [];
end

% [P12b] Constant-velocity DoA predictor. OPT-IN: with the adapt.predict.cv
% block absent this is inert and every pre-P12b result stands unchanged. With
% it present, adapt_predict_update steers the null at the CV-predicted angle
% while the jammer is drifting — see adapt_cv_init.m for the measurement that
% motivates it and for the mirror-fold problem it has to solve.
if isfield(p, 'cv') && ~isempty(p.cv)
    if nargin < 9 || isempty(sim_config)
        error('adapt_predict_init:MissingSimConfig', ...
            ['adapt.predict.cv is configured, so adapt_predict_init needs the ' ...
             'sim config (9th argument) for its dt_s. Callers that do not use ' ...
             'the CV predictor may omit it.']);
    end
    if ~isfield(sim_config, 'dt_s') || isempty(sim_config.dt_s)
        error('adapt_predict_init:MissingKey', ...
            'adapt.predict.cv needs sim.dt_s to build its motion model.');
    end
    state.cv_enabled = true;
    state.cv = adapt_cv_init(p.cv, sim_config.dt_s, E1, E2, theta_deg, phi_deg);
else
    state.cv_enabled = false;
    state.cv = [];
end

% Presence history + on/off period-detection state. The periodogram analyses
% the most recent buffer_len samples of this growing 0/1 record.
if ~state.onoff_enabled
    state.buffer_len = p.buffer_len;    % (the onoff block sizes it above)
end
state.presence    = [];                            % 1 = jammer detected, per step
state.min_periods = p.min_periods;
state.lead_steps  = p.lead_steps;
state.period_est  = NaN;                           % [steps], NaN until learned
state.duty_est    = NaN;                           % ON fraction, NaN until learned
state.last_on_k   = NaN;                           % step index of last 0->1 turn-on
state.steps_absent = 0;                            % [O] consecutive steps with no jammer
state.last_doa    = struct('theta_deg', NaN, 'phi_deg', NaN, 'idx', NaN);
state.k           = 0;

% Per-step diagnostics for closed_loop_run / plot_doa_waterfall.
state.last = struct('theta_j_deg', NaN, 'phi_j_deg', NaN, 'present', false, ...
    'predicted_on', false, 'period_est', NaN);

% Initial weights: quiescent MVDR (no jammer seen yet).
state.w = adapt_lcmv_null(state.R_hat, e_s, [], state.loading);
end
