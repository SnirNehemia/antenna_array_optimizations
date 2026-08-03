function state = adapt_freq_init(adapt_config)
% ADAPT_FREQ_INIT  Initialize the jammer carrier-frequency tracker [P10].
%
%   state = ADAPT_FREQ_INIT(adapt_config)
%
%   Companion to adapt_freq_update. Holds the smoothing / lock state around the
%   stateless per-block estimator adapt_freq_estimate — the same split as
%   adapt_predict_* around adapt_music_doa.
%
%   This tracker is a SIDE CHANNEL: it commands the RF notch centre frequency,
%   not the array weights. It therefore does not follow the
%   [w, state] = <alg>_update(state, obs) contract; it returns
%   [f_notch_hz, state] instead and runs ALONGSIDE whichever Mode C beamformer
%   is active. It consumes obs.snapshots / obs.fs_hz only, so the Mode C / Mode
%   S boundary is unaffected.
%
%   Inputs:
%       adapt_config : the config adapt section. Required sub-struct freq with
%                      keys nfft_factor, presence_snr_db, smoothing_lambda.
%                      smoothing_lambda is the weight on the RETAINED estimate
%                      (as forgetting_lambda is for R_hat): 0 = trust each block
%                      outright, ->1 = very heavy smoothing. Because a single
%                      block already lands within a few hundred Hz, smoothing
%                      here buys outlier rejection rather than accuracy.
%
%   Outputs:
%       state : struct with fields
%           cfg              : the adapt.freq sub-struct (passed to the estimator).
%           smoothing_lambda : EMA weight on the retained estimate.
%           f_hat_hz         : tracked carrier offset [Hz]; NaN until first lock.
%           have_lock        : logical, true once a line has ever been seen.
%           k                : step counter.
%           last             : diagnostics {f_hat_hz, present, snr_db, holding}
%                              refreshed every update, for closed_loop_run's log
%                              and plot_freq_waterfall.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P10].

if ~isfield(adapt_config, 'freq') || isempty(adapt_config.freq)
    error('adapt_freq_init:MissingKey', ...
        'Missing required adapt config key: ''freq'' (needs nfft_factor, presence_snr_db, smoothing_lambda).');
end
cfg = adapt_config.freq;

req_field(cfg, 'nfft_factor',      'adapt.freq');
req_field(cfg, 'presence_snr_db',  'adapt.freq');
lambda = req_field(cfg, 'smoothing_lambda', 'adapt.freq');
if lambda < 0 || lambda >= 1
    error('adapt_freq_init:BadSmoothing', ...
        'adapt.freq.smoothing_lambda must be in [0, 1) (got %g).', lambda);
end

state = struct( ...
    'cfg',              cfg, ...
    'smoothing_lambda', lambda, ...
    'f_hat_hz',         NaN, ...
    'have_lock',        false, ...
    'k',                0, ...
    'last',             struct('f_hat_hz', NaN, 'present', false, ...
                               'snr_db', NaN, 'holding', false));
end


function v = req_field(s, key, section)
% Fetch a required key or raise a descriptive error (no silent defaults).
if ~isfield(s, key) || isempty(s.(key))
    error('adapt_freq_init:MissingKey', ...
        'Missing required %s config key: ''%s''.', section, key);
end
v = s.(key);
end
