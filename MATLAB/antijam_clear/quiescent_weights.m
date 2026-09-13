function weights = quiescent_weights(signal_steering)
% ══════════════════════════════════════════════════════════════════
% QUIESCENT_WEIGHTS
% The plain beam: point at the signal, ignore everything else.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   weights = QUIESCENT_WEIGHTS(signal_steering)
%
%   The matched filter, w = e_s, scaled so that the array has exactly unit gain
%   on the signal:  w' * e_s = 1. This maximises signal-to-NOISE ratio and knows
%   nothing about any jammer.
%
%   It serves three purposes in this folder, which is why it is its own file:
%       1. the baseline -- "what you get by doing nothing";
%       2. the beam whose beamwidth defines the guard sector (array_profile);
%       3. the honest fallback when an array has no degrees of freedom left to
%          null with, or when no jammer has been detected.
%
%   Inputs:
%       signal_steering : (n_elements x 1) complex array response in the signal
%                         direction, from steering_vector. Units: V/m.
%
%   Outputs:
%       weights : (n_elements x 1) complex weights, unit gain on the signal.
%                 Units: dimensionless (1/(V/m)).

signal_steering = signal_steering(:);

% Unit-gain normalisation: w = e / (e' * e)  =>  w' * e = 1.
signal_power = real(signal_steering' * signal_steering);
if signal_power <= 0
    error('quiescent_weights:DeadDirection', ...
        ['The array radiates zero power in the signal direction, so no beam ' ...
         'can be formed there. Check the requested signal angle against the ' ...
         'element patterns.']);
end
weights = signal_steering / signal_power;
end
