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

% NaN in -> NaN out. The clamp below cannot be written as a bare
% min(max(x,-1),1): MATLAB's min/max IGNORE NaN, so max(NaN,-1) is -1 and a NaN
% input silently became acos(-1) = exactly 180 deg — a plausible-looking
% "diametrically opposite" answer. [P12] That cost a wrong finding: averaging
% this over steps where a DoA estimate was NaN (jammer not detected) read as a
% systematic antipodal lock in the MUSIC estimator that did not exist.
bad     = isnan(cos_sep);
cos_sep = min(max(cos_sep, -1), 1);   % clamp for floating-point round-off
sep_deg = rad2deg(acos(cos_sep));
sep_deg(bad) = NaN;
end
