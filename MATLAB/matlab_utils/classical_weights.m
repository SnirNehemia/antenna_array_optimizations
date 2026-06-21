function weights_complex = classical_weights(name, n_side, d_over_lambda, ...
        steer_theta_deg, steer_phi_deg, element_patterns_complex, theta_deg, phi_deg)
% CLASSICAL_WEIGHTS  Complex weights for a classical tapering + steering method.
%
%   Separable outer-product aperture taper (two identical 1-D windows) times a
%   progressive steering phase. Steering mode is data-driven when element
%   patterns + grids are supplied, otherwise the geometric model is used.
%
%   Inputs:
%       name                          : classical technique (see window_1d).
%       n_side, d_over_lambda         : array geometry.
%       steer_theta_deg, steer_phi_deg: steering direction. Units: degrees.
%       element_patterns_complex      : optional (N_elements,N_theta,N_phi);
%                                       triggers data-driven steering. Default [].
%       theta_deg, phi_deg            : optional grids for data-driven steering.
%
%   Outputs:
%       weights_complex : (n_side^2 x 1) complex.
%
%   Ported from scripts/compare_classical.py:classical_weights.
%   Part of: Antenna Array Pattern Optimization Tool.

if nargin < 6, element_patterns_complex = []; end
if nargin < 7, theta_deg = []; end
if nargin < 8, phi_deg   = []; end

window_1d_vec = window_1d(name, n_side);                 % (n_side x 1)
window_2d     = window_1d_vec * window_1d_vec.';          % (n_side x n_side)

% Flatten in C-order (row-major) to match the build_ura element ordering.
% [PORT] Python window_2d.ravel(); symmetric here so order is immaterial.
weights_taper_complex = complex(reshape(window_2d.', [], 1));   % (n_side^2 x 1)

use_data_driven = ~isempty(element_patterns_complex) ...
    && ~isempty(theta_deg) && ~isempty(phi_deg);
if use_data_driven
    steering_phases = data_driven_steering_vector( ...
        element_patterns_complex, theta_deg, phi_deg, steer_theta_deg, steer_phi_deg);
else
    steering_phases = steering_phase_vector( ...
        n_side, d_over_lambda, steer_theta_deg, steer_phi_deg);
end

weights_complex = weights_taper_complex .* steering_phases;
end
