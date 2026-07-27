function [it, ip] = nearest_index_2d(theta_deg, phi_deg, target_theta_deg, target_phi_deg)
% NEAREST_INDEX_2D  Grid indices closest to a target (theta, phi) point.
%
%   [it, ip] = NEAREST_INDEX_2D(theta_deg, phi_deg, target_theta_deg, target_phi_deg)
%
%   The far-field grid is rectangular (independent theta/phi axes), so the
%   nearest grid point is exactly the pair of per-axis nearest indices.
%
%   Inputs:
%       theta_deg, phi_deg : grid axis vectors. Units: degrees.
%       target_theta_deg, target_phi_deg : target point (scalars). Units: degrees.
%
%   Outputs:
%       it, ip : 1-based indices into theta_deg / phi_deg (MATLAB convention).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P7].

it = nearest_index(theta_deg(:), target_theta_deg);
ip = nearest_index(phi_deg(:), mod(target_phi_deg, 360));
end
