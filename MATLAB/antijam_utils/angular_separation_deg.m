function sep_deg = angular_separation_deg(theta1_deg, phi1_deg, theta2_deg, phi2_deg)
% ANGULAR_SEPARATION_DEG  True spherical angular distance between two (theta, phi) points.
%
%   sep_deg = ANGULAR_SEPARATION_DEG(theta1_deg, phi1_deg, theta2_deg, phi2_deg)
%
%   theta is measured from boresight (0 = broadside), phi is azimuth around
%   boresight — the same (theta, phi) convention as the CST far-field grid.
%   Inputs may be scalars or same-shape arrays (elementwise).
%
%   Inputs:
%       theta1_deg, phi1_deg : first point(s). Units: degrees.
%       theta2_deg, phi2_deg : second point(s). Units: degrees.
%
%   Outputs:
%       sep_deg : great-circle angular separation, in [0, 180]. Units: degrees.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P7].

t1 = deg2rad(theta1_deg);
t2 = deg2rad(theta2_deg);
dphi = deg2rad(phi1_deg - phi2_deg);

cos_sep = cos(t1) .* cos(t2) + sin(t1) .* sin(t2) .* cos(dphi);
cos_sep = min(max(cos_sep, -1), 1);   % clamp for floating-point round-off
sep_deg = rad2deg(acos(cos_sep));
end
