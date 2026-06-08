function cut_vm = principal_plane_cut(array_factor_complex, phi_deg)
% PRINCIPAL_PLANE_CUT  Principal-plane cut from a 2-D array factor.
%
%   Combines the phi=0 half-plane (positive-theta side) and phi=180 half-plane
%   (negative-theta side) into one cut spanning -90..+90 degrees:
%       theta_display < 0  -> phi=180 half-plane
%       theta_display > 0  -> phi=0   half-plane
%
%   Inputs:
%       array_factor_complex : (N_theta x N_phi) complex [V/m] (or real magnitude).
%       phi_deg              : azimuth grid (N_phi). Units: degrees.
%
%   Outputs:
%       cut_vm : (2*N_theta-1 x 1) linear magnitude [V/m].
%
%   Ported from scripts/compare_classical.py:principal_plane_cut.
%   Part of: Antenna Array Pattern Optimization Tool.

phi0_idx   = nearest_index(phi_deg, 0.0);
phi180_idx = nearest_index(phi_deg, 180.0);

af_pos_vm = abs(array_factor_complex(:, phi0_idx));     % phi=0
af_neg_vm = abs(array_factor_complex(:, phi180_idx));   % phi=180

% [PORT] np.concatenate([af_neg[::-1][:-1], af_pos]).
neg_rev = flipud(af_neg_vm);
cut_vm  = [neg_rev(1:end - 1); af_pos_vm];
end
