function theta_display_deg = principal_plane_theta_axis(theta_deg)
% PRINCIPAL_PLANE_THETA_AXIS  Symmetric theta display axis for a principal cut.
%
%   Mirrors the one-sided theta grid [0,90] into [-90,...,0,...,+90] of length
%   2*N_theta-1.
%
%   Inputs:
%       theta_deg : one-sided elevation grid (N_theta), starting at 0. Units: deg.
%
%   Outputs:
%       theta_display_deg : (2*N_theta-1 x 1). Units: degrees.
%
%   Ported from scripts/compare_classical.py:principal_plane_theta_axis.
%   Part of: Antenna Array Pattern Optimization Tool.

% [PORT] np.concatenate([-theta_deg[::-1][:-1], theta_deg]).
rev = flipud(theta_deg(:));
theta_display_deg = [-rev(1:end - 1); theta_deg(:)];
end
