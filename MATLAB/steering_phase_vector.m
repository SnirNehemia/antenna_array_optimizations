function steering_phases = steering_phase_vector(n_side, d_over_lambda, steer_theta_deg, steer_phi_deg)
% STEERING_PHASE_VECTOR  Geometric 2-D steering phase vector for a URA.
%
%   Conjugate of the array manifold at the steering direction; multiplying
%   element weights by this vector points the beam to (steer_theta, steer_phi).
%
%   Inputs:
%       n_side, d_over_lambda             : array geometry.
%       steer_theta_deg, steer_phi_deg    : steering direction. Units: degrees.
%
%   Outputs:
%       steering_phases : (n_side^2 x 1) complex, unit amplitude.
%
%   Ported from scripts/compare_classical.py:_steering_phase_vector.
%   Part of: Antenna Array Pattern Optimization Tool.

steer_theta_rad = deg2rad(steer_theta_deg);
steer_phi_rad   = deg2rad(steer_phi_deg);
kd              = 2 * pi * d_over_lambda;
offset          = (n_side - 1) / 2;

u0 = sin(steer_theta_rad) * cos(steer_phi_rad);
v0 = sin(steer_theta_rad) * sin(steer_phi_rad);

steering_phases = complex(zeros(n_side * n_side, 1));
elem_idx = 0;
for m = 0:(n_side - 1)
    for nn = 0:(n_side - 1)
        m_c = m  - offset;
        n_c = nn - offset;
        phase_rad = kd * (m_c * u0 + n_c * v0);
        elem_idx  = elem_idx + 1;
        steering_phases(elem_idx) = exp(-1i * phase_rad);
    end
end
end
