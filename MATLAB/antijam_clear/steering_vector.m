function steering = steering_vector(element_patterns, theta_deg, phi_deg, ...
                                    look_theta_deg, look_phi_deg)
% ══════════════════════════════════════════════════════════════════
% STEERING_VECTOR
% The array's complex response to a plane wave from one direction.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   steering = STEERING_VECTOR(element_patterns, theta_deg, phi_deg, ...
%                              look_theta_deg, look_phi_deg)
%
%   THE physical primitive of this folder. A plane wave of complex amplitude a
%   arriving from (theta, phi) puts
%
%       x = a * e(theta, phi)
%
%   across the elements, where e is this function's output. Everything else is
%   written in terms of e: it is the direction you point at (signal), the
%   direction you cancel (jammer), the axis MUSIC scans, and the vector the
%   oracle is handed. One definition used four ways.
%
%   CONJUGATION CONVENTION -- read this once and it holds everywhere.
%   The array output is
%
%       y = w' * x            (w' is the Hermitian transpose)
%
%   so the array's complex gain in direction e is  w' * e, and the output power
%   in that direction is |w' * e|^2. Every formula in this folder uses that form
%   without exception.
%
%   matlab_utils/compute_array_factor uses the OTHER convention -- it forms
%   sum_n w_n * E_n, which is w.' * e, with no conjugate. The two differ by
%   conj(w), and that conjugate is applied at exactly the two places this folder
%   calls compute_array_factor, both of them commented:
%
%       array_profile    -- measuring the quiescent beam's shape;
%       plot_pattern_cut -- drawing the pattern.
%
%   Nowhere else. Keep it that way: a beamformer that is right in the table and
%   mirrored in the figure is very expensive to diagnose.
%
%   (Two other uses of conj() in this folder are unrelated to the convention:
%   the elementwise form of e' * R * e in estimate_jammer_angle, and a scalar
%   normalisation in max_sinr_weights / oracle_weights.)
%
%   NO NORMALISATION. The returned vector carries the element patterns' true
%   relative magnitudes, so a direction where the elements radiate weakly gives
%   a short vector -- which is physically correct and is what makes the
%   achievable SINR array-dependent. (Contrast matlab_utils/
%   data_driven_steering_vector, which returns conj(e)./|e| -- a unit-modulus
%   matched-filter WEIGHT, not a steering vector. Different object, different
%   convention; it is deliberately not used here.)
%
%   NEAREST GRID POINT. The look direction is snapped to the CST grid (1 deg
%   here) rather than interpolated. The resulting angular error is bounded by
%   half a grid step and is a real, honest contributor to steering mismatch --
%   the effect diagonal loading exists to absorb.
%
%   [FUTURE] Two-component operation: this returns (n_elements x 2), one column
%   per polarization, and every |w' * e|^2 downstream becomes
%   sum(abs(w' * e).^2) / n_components.
%
%   Inputs:
%       element_patterns : (n_elements x n_theta x n_phi) complex. Units: V/m.
%       theta_deg        : (n_theta x 1) elevation grid. Units: degrees.
%       phi_deg          : (n_phi x 1)   azimuth grid.   Units: degrees.
%       look_theta_deg   : look-direction elevation. Units: degrees.
%       look_phi_deg     : look-direction azimuth.   Units: degrees.
%
%   Outputs:
%       steering : (n_elements x 1) complex array response. Units: V/m.

% ────────────────────────── SNAP TO GRID ──────────────────────────

theta_index = nearest_index(theta_deg, look_theta_deg);
phi_index   = nearest_index(phi_deg,   look_phi_deg);

% ────────────────────────── READ THE MANIFOLD ─────────────────────

% One column of the pattern cube: every element's field in this direction.
steering = element_patterns(:, theta_index, phi_index);
steering = steering(:);
end
