function [f_notch_hz, state] = adapt_freq_update(state, obs)
% ADAPT_FREQ_UPDATE  Track the jammer carrier; command the RF notch [P10].
%
%   [f_notch_hz, state] = ADAPT_FREQ_UPDATE(state, obs)
%
%   One step of the carrier tracker. Runs adapt_freq_estimate on this step's
%   snapshot block, then folds the result into the tracked estimate with three
%   behaviours:
%
%     LOCK    — the first block with a detected line sets the estimate outright.
%               Smoothing from an initial NaN would never converge.
%     SMOOTH  — subsequent detections are EMA-blended,
%               f <- lambda*f + (1-lambda)*f_new, rejecting the occasional
%               outlier block without materially slowing the tracker.
%     HOP     — a detection further than HOP_RESET_BINS periodogram bins from
%               the tracked value is treated as a genuine carrier change, not
%               estimator jitter, and REPLACES the estimate instead of being
%               averaged with it. Without this, smoothing would walk the notch
%               slowly across the gap after a hop, parking it on neither
%               carrier for many steps. The threshold is derived from the FFT
%               bin width rather than configured: sub-bin jitter is a small
%               fraction of one bin, so four bins is unambiguous either way.
%
%   HOLD-THROUGH-OFF: when no line is detected — the jammer is silent, or the
%   block was unlucky — the tracked estimate is HELD, not cleared. The notch
%   therefore stays parked on the last known carrier through an OFF gap and is
%   already correct the instant the jammer returns. This is the spectral
%   analogue of adapt_predict_update holding last_doa through OFF gaps, and it
%   is why the notch costs its (tiny) insertion loss continuously rather than
%   only while the jammer transmits.
%
%   MODE C CONTRACT: consumes obs.snapshots and obs.fs_hz only. Never touches
%   obs.sinr_db, and has no access to the true carrier.
%
%   Inputs:
%       state : struct from adapt_freq_init (advanced state returned).
%       obs   : observation struct from sim_engine_step. Requires non-empty
%               snapshots and fs_hz (i.e. Mode C with the waveform layer on).
%
%   Outputs:
%       f_notch_hz : commanded notch centre, as an offset from the array centre
%                    frequency [Hz]. NaN before the first lock — sim_engine_step
%                    reads NaN as "notch not engaged", so the notch costs
%                    nothing until the tracker actually knows something.
%       state      : updated state.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P10].

HOP_RESET_BINS = 4;   % see HOP above; well outside sub-bin estimator jitter

if ~isfield(obs, 'snapshots') || isempty(obs.snapshots)
    error('adapt_freq_update:NoSnapshots', ...
        'adapt_freq_update is a Mode C consumer and needs obs.snapshots.');
end
if ~isfield(obs, 'fs_hz') || isempty(obs.fs_hz)
    error('adapt_freq_update:NoSampleRate', ...
        ['obs.fs_hz is empty — the [P10] waveform layer is off, so the snapshots ' ...
         'are spectrally white and carry no carrier to estimate. Set sim.fs_hz.']);
end

state.k = state.k + 1;
freq    = adapt_freq_estimate(obs.snapshots, obs.fs_hz, state.cfg);

holding = true;
if freq.present
    holding = false;
    df = obs.fs_hz / numel(freq.f_axis_hz);          % periodogram bin width
    if ~state.have_lock
        state.f_hat_hz  = freq.f_hat_hz;             % LOCK
        state.have_lock = true;
    elseif abs(freq.f_hat_hz - state.f_hat_hz) > HOP_RESET_BINS * df
        state.f_hat_hz = freq.f_hat_hz;              % HOP
    else
        lam = state.smoothing_lambda;                % SMOOTH
        state.f_hat_hz = lam * state.f_hat_hz + (1 - lam) * freq.f_hat_hz;
    end
end

% pspec / f_axis_hz are OVERWRITTEN each step, never accumulated — callers that
% want a history (plot_freq_waterfall) sample them on their own stride.
state.last = struct('f_hat_hz', state.f_hat_hz, 'present', freq.present, ...
                    'snr_db', freq.snr_db, 'holding', holding, ...
                    'pspec', freq.pspec, 'f_axis_hz', freq.f_axis_hz);
f_notch_hz = state.f_hat_hz;
end
