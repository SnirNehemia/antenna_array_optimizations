function weights_norm = power_normalize_weights(weights_complex)
% POWER_NORMALIZE_WEIGHTS  Normalise complex weights to unit total power.
%
%   Divides by the L2 norm of the weight vector (matches the internal
%   normalisation in build_cost_function), so a fully coherent uniform array
%   peaks at sqrt(N) -> 10*log10(N) dBi.
%
%   Inputs:
%       weights_complex : (N x 1) complex. Units: dimensionless.
%
%   Outputs:
%       weights_norm : (N x 1) complex, power-normalised.
%
%   Ported from scripts/compare_classical.py:_power_normalize_weights.
%   Part of: Antenna Array Pattern Optimization Tool.

total_power  = sum(abs(weights_complex(:)) .^ 2);
safe_power   = max(total_power, 1e-30);
weights_norm = weights_complex(:) / sqrt(safe_power);
end
