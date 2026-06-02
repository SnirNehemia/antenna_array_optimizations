function directivity_dbi_grid = compute_directivity_dbi_grid(array_factor, theta_deg, phi_deg, normalizer_power)
% COMPUTE_DIRECTIVITY_DBI_GRID  Directivity in dBi over the full angular grid.
%
%   directivity_dbi_grid = COMPUTE_DIRECTIVITY_DBI_GRID(array_factor, theta_deg, ...
%                                                       phi_deg, normalizer_power)
%
%   Discrete rectangle-rule spherical integral:
%       P_total ~ sum |AF|^2 sin(theta) dtheta dphi
%       D(theta,phi) = 4*pi*|AF|^2 / P_total
%
%   Inputs:
%       array_factor     : complex AF, (N_theta x N_phi). Units: V/m (arbitrary scale).
%       theta_deg        : elevation grid, (N_theta x 1). Units: degrees (0..180).
%       phi_deg          : azimuth grid,   (N_phi x 1).   Units: degrees.
%       normalizer_power : optional precomputed total radiated power as the
%                          denominator (e.g. P_copol + P_cross). Pass [] to use
%                          the power radiated by array_factor alone. Default [].
%
%   Outputs:
%       directivity_dbi_grid : (N_theta x N_phi). Units: dBi.
%
%   Ported from src/metrics/metrics.py:_compute_directivity_dbi_grid.
%   Part of: Antenna Array Pattern Optimization Tool.

LOG10_EPSILON = 1e-30;
if nargin < 4
    normalizer_power = [];
end

power_grid = abs(array_factor) .^ 2;

theta_rad = deg2rad(theta_deg(:));
if numel(theta_rad) > 1
    dtheta_rad = mean(diff(theta_rad));
else
    dtheta_rad = pi;
end
if numel(phi_deg) > 1
    dphi_rad = deg2rad(mean(diff(phi_deg(:))));
else
    dphi_rad = 2 * pi;
end

sin_theta = sin(theta_rad);   % (N_theta x 1), solid-angle weighting

if ~isempty(normalizer_power)
    p_total = double(normalizer_power);
else
    % Numerical spherical integral: sum |AF|^2 sin(theta) dtheta dphi.
    p_total = sum(power_grid .* sin_theta, 'all') * dtheta_rad * dphi_rad;
end

directivity_linear = (4 * pi * power_grid) / max(p_total, 1e-30);
directivity_dbi_grid = 10 * log10(directivity_linear + LOG10_EPSILON);
end
