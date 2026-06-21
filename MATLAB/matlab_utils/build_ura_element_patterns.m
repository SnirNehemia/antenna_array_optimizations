function element_patterns_complex = build_ura_element_patterns(n_side, d_over_lambda, theta_deg, phi_deg)
% BUILD_URA_ELEMENT_PATTERNS  Synthetic isotropic N x N URA element-pattern stack.
%
%   ep = BUILD_URA_ELEMENT_PATTERNS(n_side, d_over_lambda, theta_deg, phi_deg)
%
%   Each element is an isotropic radiator with the free-space progressive phase
%       E_{mn}(theta,phi) = exp(1i k d (m_c sin(theta) cos(phi) + n_c sin(theta) sin(phi)))
%   with centre-referenced indices m_c, n_c.
%
%   Inputs:
%       n_side        : elements per side; total = n_side^2.
%       d_over_lambda : element spacing in wavelengths. Units: dimensionless.
%       theta_deg     : elevation grid (N_theta). Units: degrees.
%       phi_deg       : azimuth grid (N_phi).     Units: degrees.
%
%   Outputs:
%       element_patterns_complex : (n_side^2, N_theta, N_phi) complex [V/m].
%
%   Ported from scripts/compare_classical.py:build_ura_element_patterns.
%   Part of: Antenna Array Pattern Optimization Tool.

n_elements = n_side * n_side;
n_theta    = numel(theta_deg);
n_phi      = numel(phi_deg);

theta_rad = deg2rad(theta_deg(:));     % (N_theta x 1)
phi_rad   = deg2rad(phi_deg(:));       % (N_phi x 1)

sin_theta = sin(theta_rad);
cos_phi   = cos(phi_rad);
sin_phi   = sin(phi_rad);
u_grid    = sin_theta * cos_phi.';     % (N_theta x N_phi)
v_grid    = sin_theta * sin_phi.';

kd     = 2 * pi * d_over_lambda;
offset = (n_side - 1) / 2;

element_patterns_complex = complex(zeros(n_elements, n_theta, n_phi));
elem_idx = 0;
for m = 0:(n_side - 1)
    for nn = 0:(n_side - 1)
        m_c = m  - offset;
        n_c = nn - offset;
        phase_rad = kd * (m_c * u_grid + n_c * v_grid);   % (N_theta x N_phi)
        elem_idx  = elem_idx + 1;
        element_patterns_complex(elem_idx, :, :) = exp(1i * phase_rad);
    end
end
end
