function [h2, loss_frac] = sim_notch_response(notch_cfg, f_notch_hz, f_probe_hz, fs_hz)
% SIM_NOTCH_RESPONSE  RF notch power response + wideband insertion loss [P10].
%
%   [h2, loss_frac] = SIM_NOTCH_RESPONSE(notch_cfg, f_notch_hz, f_probe_hz, fs_hz)
%
%   Closed-form second-order notch, centred at f_notch_hz:
%
%       d       = 10^(-depth_db/10)                       % floor at the centre
%       a       = bw_hz / 2                               % half 3-dB width
%       |H(f)|^2 = (D^2 + d*a^2) / (D^2 + a^2),  D = f - f_notch_hz
%
%   so |H|^2 = d exactly at the centre and (1+d)/2 at f = f_notch +- bw_hz/2 —
%   i.e. bw_hz really is the 3-dB width.
%
%   WHY A SCALAR MODEL IS EXACT HERE: the notch is a linear time-invariant
%   filter applied identically on every element, so it COMMUTES with the linear
%   beamformer. Filtering each snapshot sample and then combining gives exactly
%   the same output powers as scaling the combined powers, which is all
%   sim_engine_step's pattern-level SINR needs. No per-sample filtering, no
%   filter coefficients, no group-delay bookkeeping.
%
%   The two outputs enter the SINR equation on opposite sides:
%     * h2        multiplies the JAMMER power. It is evaluated at the jammer's
%                 TRUE carrier while the notch sits at the ESTIMATED one, so a
%                 frequency-estimation error directly costs rejection — this is
%                 what couples the P10 accuracy KPI to delivered SINR.
%     * loss_frac multiplies the DESIRED SIGNAL and the NOISE, which are both
%                 white across the captured band (the confirmed wideband-signal
%                 assumption), so each loses the same integrated fraction.
%
%   loss_frac uses the exact band-limited integral of 1 - |H|^2 over
%   [-fs_hz/2, fs_hz/2] relative to f_notch_hz,
%
%       loss_frac = (1-d) * a * [atan((fs/2 - f_n)/a) - atan((-fs/2 - f_n)/a)] / fs
%
%   which reduces to the familiar pi*(1-d)*bw_hz / (2*fs_hz) whenever the notch
%   sits well inside the band. At bw_hz = 200 kHz, fs_hz = 24 MHz, depth 35 dB
%   that is 1.3% of the signal power — 0.057 dB — bought for 35 dB of jammer
%   rejection. The notch is very nearly free against a wideband desired signal.
%
%   Inputs:
%       notch_cfg  : struct. Required: depth_db (>= 0), bw_hz (> 0).
%                    Field 'mode' is interpreted by the CALLER (sim_engine_step
%                    picks f_notch_hz from it); this function is pure response.
%       f_notch_hz : notch centre, as an offset from the array centre frequency
%                    [Hz] (baseband convention, same as f_probe_hz).
%       f_probe_hz : frequency at which to evaluate |H|^2 [Hz]. May be an array.
%       fs_hz      : captured bandwidth / snapshot sample rate [Hz].
%
%   Outputs:
%       h2        : |H(f_probe_hz)|^2, linear power gain in [d, 1]. Same size as
%                   f_probe_hz.
%       loss_frac : scalar in [0, 1], fraction of a band-filling white signal's
%                   power removed by the notch.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P10].

depth_db = req_field(notch_cfg, 'depth_db', 'notch');
bw_hz    = req_field(notch_cfg, 'bw_hz',    'notch');

if depth_db < 0
    error('sim_notch_response:BadDepth', ...
        'notch.depth_db must be >= 0 (got %g); it is an attenuation in dB.', depth_db);
end
if bw_hz <= 0
    error('sim_notch_response:BadBandwidth', ...
        'notch.bw_hz must be > 0 (got %g).', bw_hz);
end
if fs_hz <= 0
    error('sim_notch_response:BadSampleRate', 'fs_hz must be > 0 (got %g).', fs_hz);
end

d = 10^(-depth_db / 10);
a = bw_hz / 2;

delta = f_probe_hz - f_notch_hz;
h2    = (delta.^2 + d * a^2) ./ (delta.^2 + a^2);

% Exact band-limited integral of the stopped fraction (1 - |H|^2). Using the
% finite band rather than the infinite-line approximation keeps loss_frac <= 1
% even for a pathologically wide notch.
hi = ( fs_hz / 2 - f_notch_hz) / a;
lo = (-fs_hz / 2 - f_notch_hz) / a;
loss_frac = (1 - d) * a * (atan(hi) - atan(lo)) / fs_hz;
loss_frac = min(max(loss_frac, 0), 1);
end


function v = req_field(s, key, section)
% Fetch a required key or raise a descriptive error (no silent defaults).
if ~isfield(s, key) || isempty(s.(key))
    error('sim_notch_response:MissingKey', ...
        'Missing required %s config key: ''%s''.', section, key);
end
v = s.(key);
end
