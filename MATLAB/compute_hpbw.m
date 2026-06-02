function width = compute_hpbw(cut_dbi, peak_idx, opposite_cut)
% COMPUTE_HPBW  3 dB beamwidth (in grid-index units) by outward scan from a peak.
%
%   width = COMPUTE_HPBW(cut_dbi, peak_idx, opposite_cut)
%
%   Scans outward from peak_idx with linear interpolation to the -3 dB crossing.
%   When a scan reaches a grid boundary without crossing and opposite_cut (the
%   cut at azimuth phi+180) is supplied, the beam is assumed to wrap around the
%   pole and the search continues in that cut.
%
%   Inputs:
%       cut_dbi      : 1-D directivity cut, dBi.
%       peak_idx     : 0-BASED index of the global peak within cut_dbi (mirrors
%                      the Python convention; callers pass MATLAB index minus 1).
%       opposite_cut : optional theta-cut at phi+180 for pole wrap-around. Pass
%                      [] to disable (e.g. for phi-cuts). Default [].
%
%   Outputs:
%       width : full 3 dB beamwidth in grid-index units (multiply by the angular
%               step outside to convert to degrees).
%
%   Ported from src/metrics/metrics.py:_compute_hpbw.
%   Part of: Antenna Array Pattern Optimization Tool.

if nargin < 3
    opposite_cut = [];
end
cut_dbi   = cut_dbi(:);
n         = numel(cut_dbi);
threshold = cut_dbi(peak_idx + 1) - 3.0;

% ── Left (decreasing-index) crossing ──────────────────────────────
left       = peak_idx;
left_found = false;
for i = peak_idx - 1 : -1 : 0
    if cut_dbi(i + 1) < threshold
        denom = cut_dbi(i + 2) - cut_dbi(i + 1);          % > 0
        if abs(denom) <= 1e-30, denom = 1e-30; end
        t    = (threshold - cut_dbi(i + 1)) / denom;
        left = i + clip01(t);
        left_found = true;
        break;
    end
end
if ~left_found && ~isempty(opposite_cut)
    oc = opposite_cut(:);
    for i = 1 : n - 1
        if oc(i + 1) < threshold
            denom = oc(i + 1) - oc(i);                     % < 0
            if abs(denom) <= 1e-30, denom = -1e-30; end
            t         = (threshold - oc(i)) / denom;
            opp_extra = (i - 1) + clip01(t);
            left      = -opp_extra;
            break;
        end
    end
end

% ── Right (increasing-index) crossing ─────────────────────────────
right       = peak_idx;
right_found = false;
for i = peak_idx + 1 : n - 1
    if cut_dbi(i + 1) < threshold
        denom = cut_dbi(i + 1) - cut_dbi(i);              % < 0
        if abs(denom) <= 1e-30, denom = -1e-30; end
        t     = (threshold - cut_dbi(i)) / denom;
        right = (i - 1) + clip01(t);
        right_found = true;
        break;
    end
end
if ~right_found && ~isempty(opposite_cut)
    oc = opposite_cut(:);
    for i = n - 2 : -1 : 0
        if oc(i + 1) < threshold
            denom = oc(i + 2) - oc(i + 1);                 % > 0
            if abs(denom) <= 1e-30, denom = 1e-30; end
            t         = (threshold - oc(i + 1)) / denom;
            opp_extra = (n - 1) - (i + clip01(t));
            right     = (n - 1) + opp_extra;
            break;
        end
    end
end

width = right - left;
end


function y = clip01(t)
% Clip to [0, 1] (np.clip(t, 0.0, 1.0)).
y = min(max(t, 0.0), 1.0);
end
