function steering = data_driven_steering_vector(element_patterns_complex, theta_deg, phi_deg, ...
                                                steer_theta_deg, steer_phi_deg)
% DATA_DRIVEN_STEERING_VECTOR  Model-free steering vector from element-pattern data.
%
%   Evaluates each element's far-field at the nearest grid point to the target
%   direction and returns the conjugate unit-phase vector. Correct for real CST
%   patterns; identical to steering_phase_vector for synthetic phase-only patterns.
%
%   Inputs:
%       element_patterns_complex      : (N_elements, N_theta, N_phi) complex [V/m].
%       theta_deg, phi_deg            : angle grids. Units: degrees.
%       steer_theta_deg, steer_phi_deg: steering direction. Units: degrees.
%
%   Outputs:
%       steering : (N_elements x 1) unit-amplitude complex.
%
%   Ported from scripts/compare_classical.py:_data_driven_steering_vector.
%   Part of: Antenna Array Pattern Optimization Tool.

theta_idx = nearest_index(theta_deg, steer_theta_deg);
phi_idx   = nearest_index(phi_deg,   steer_phi_deg);

% Array manifold: element responses at the steering direction.
manifold = element_patterns_complex(:, theta_idx, phi_idx);   % (N_elements x 1)
manifold = manifold(:);

amp      = max(abs(manifold), 1e-30);
steering = conj(manifold) ./ amp;
end
