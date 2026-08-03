function freq = adapt_freq_estimate(X, fs_hz, cfg)
% ADAPT_FREQ_ESTIMATE  Jammer carrier estimate from one snapshot block [P10].
%
%   freq = ADAPT_FREQ_ESTIMATE(X, fs_hz, cfg)
%
%   Pure estimator (no state), the spectral counterpart of adapt_music_doa:
%   given one Mode C snapshot block it returns the dominant spectral line's
%   frequency, a presence flag, and the periodogram for the waterfall plot.
%   State — smoothing, hold-through-OFF, hop detection — lives in
%   adapt_freq_init / adapt_freq_update, exactly as last_doa lives in
%   adapt_predict_* rather than in adapt_music_doa.
%
%   Method, and why:
%     1. INCOHERENT AVERAGE of the per-element periodograms, NOT the beamformer
%        output. Once the spatial null has formed, w'X has the jammer suppressed
%        by 30-70 dB and its carrier is unmeasurable — the estimator would go
%        blind exactly when it is working. Reading the raw element snapshots
%        keeps the frequency observable for the whole run, and needs no DoA
%        estimate (so it does not inherit MUSIC's failure modes either).
%     2. NO WINDOW (rectangular). Tapering is the reflex here and it is WRONG
%        for this signal environment: a window suppresses spectral leakage from
%        strong NARROWBAND components, and under the locked single-jammer scope
%        the only narrowband component in the band is the jammer itself. The
%        desired signal is white, so it has no sidelobes to smear. A Hann taper
%        therefore buys nothing and costs variance — measured back-to-back on
%        the toy fixture (200 trials, K = 128, fs = 24 MHz, J/N = 20 dB):
%            sigma_s_db =  3   Hann 344 Hz RMSE  vs  rectangular 242 Hz
%            sigma_s_db = 20   Hann 1247 Hz      vs  rectangular 614 Hz
%        Revisit this if a second narrowband emitter is ever added to the model;
%        until then the rectangular window is strictly better.
%     3. ZERO-PADDED FFT (nfft_factor x K) then PARABOLIC INTERPOLATION IN dB of
%        the peak. This is the same sub-bin refinement detect_period uses in
%        adapt_predict_update; the raw bin alone is far too coarse (fs/K =
%        187 kHz at fs = 24 MHz, K = 128). Zero-padding to 4x samples the main
%        lobe densely enough that a parabola fits it well near the peak; 8x was
%        measured to add nothing (343 vs 344 Hz). Note that Jacobsen-type
%        interpolators do NOT apply here — they assume adjacent un-padded bins,
%        and on a 4x-padded spectrum they degrade to ~22 kHz.
%     4. PRESENCE from the peak-to-MEDIAN ratio of the periodogram, mirroring
%        adapt_music_doa's eigengap test. Median, not mean, so the peak itself
%        does not inflate its own background reference.
%
%   ACCURACY, measured (not the textbook bound): 242 Hz RMSE at sigma_s_db = 3
%   and 614 Hz at sigma_s_db = 20, K = 128, fs = 24 MHz, J/N = 20 dB. What sets
%   it is the in-band ratio sigma_j^2|e_j|^2 / (sigma_s^2|e_s|^2 + sigma_n^2) —
%   the WIDEBAND DESIRED SIGNAL, not the noise floor, is the dominant competitor
%   once sigma_s_db approaches jn_ratio_db, which is why accuracy halves between
%   those two operating points. Either way this is well under 1% of a typical
%   200 kHz notch width, so estimation error costs about 1 dB of the notch's
%   35 dB rejection. The estimator is not the bottleneck in this phase.
%
%   Inputs:
%       X     : (N_el x K) complex snapshot block (obs.snapshots).
%       fs_hz : snapshot sample rate [Hz] (obs.fs_hz).
%       cfg   : struct. Required: nfft_factor (>= 1), presence_snr_db.
%
%   Outputs:
%       freq : struct with fields
%           present   : logical, true when a line clears presence_snr_db.
%           f_hat_hz  : estimated carrier OFFSET from the array centre
%                       frequency [Hz]; NaN when not present.
%           snr_db    : peak-to-median periodogram ratio [dB].
%           pspec     : (1 x Nfft) periodogram, linear power, fftshifted to
%                       match f_axis_hz.
%           f_axis_hz : (1 x Nfft) frequency axis [Hz], -fs/2 .. +fs/2.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P10].

nfft_factor    = req_field(cfg, 'nfft_factor',    'adapt.freq');
presence_snr_db = req_field(cfg, 'presence_snr_db', 'adapt.freq');

if isempty(X)
    error('adapt_freq_estimate:NoSnapshots', ...
        'X is empty — adapt_freq_estimate needs Mode C snapshots.');
end
K = size(X, 2);
if K < 8
    error('adapt_freq_estimate:TooFewSnapshots', ...
        'Need K >= 8 snapshots for a periodogram (got %d).', K);
end
if nfft_factor < 1
    error('adapt_freq_estimate:BadNfftFactor', ...
        'adapt.freq.nfft_factor must be >= 1 (got %g).', nfft_factor);
end
if fs_hz <= 0
    error('adapt_freq_estimate:BadSampleRate', 'fs_hz must be > 0 (got %g).', fs_hz);
end

% Rectangular window (see note 2 above — deliberate, not an oversight).
nfft = round(nfft_factor * K);

% Incoherent average of the per-element periodograms.
P = mean(abs(fft(X, nfft, 2)).^2, 1);                      % (1 x nfft)

% Centre the axis on DC so a negative carrier offset reads as negative.
P      = fftshift(P);
df     = fs_hz / nfft;
f_axis = (-floor(nfft / 2):(nfft - 1 - floor(nfft / 2))) * df;

[pk, i_pk] = max(P);
bg     = median(P);
snr_db = 10 * log10(pk / max(bg, realmin));

freq = struct('present', false, 'f_hat_hz', NaN, 'snr_db', snr_db, ...
              'pspec', P, 'f_axis_hz', f_axis);
if snr_db < presence_snr_db
    return
end

% Sub-bin refinement: parabolic vertex through the peak and its neighbours, in
% dB. Bins at the very edge of the axis cannot be refined (delta stays 0).
delta = 0;
if i_pk > 1 && i_pk < nfft
    a = 10 * log10(max(P(i_pk - 1), realmin));
    b = 10 * log10(max(P(i_pk),     realmin));
    c = 10 * log10(max(P(i_pk + 1), realmin));
    denom = a - 2 * b + c;
    if denom < 0                          % a genuine maximum, not a flat/edge
        delta = 0.5 * (a - c) / denom;
        delta = min(max(delta, -0.5), 0.5);
    end
end

freq.present  = true;
freq.f_hat_hz = f_axis(i_pk) + delta * df;
end


function v = req_field(s, key, section)
% Fetch a required key or raise a descriptive error (no silent defaults).
if ~isfield(s, key) || isempty(s.(key))
    error('adapt_freq_estimate:MissingKey', ...
        'Missing required %s config key: ''%s''.', section, key);
end
v = s.(key);
end
