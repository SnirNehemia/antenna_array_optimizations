function [theta_ext_deg, phi_ext_deg, theta_phys_idx, phi_offset, sin_theta_ext] = ...
        extended_grid_maps(theta_deg, phi_deg)
% EXTENDED_GRID_MAPS  Build the mirrored/tiled extended (theta, phi) grid and the
% maps back to physical grid indices. Shared by build_cost_function and
% build_directive_physical_masks (the Python code inlines this in both).
%
%   Theta is mirrored at both poles; phi is tiled three times so rectangular
%   windows near the poles or straddling the phi=0/360 seam are handled without
%   special cases.
%
%   Extended layout (3*N_theta-2 rows):
%       top mirror   : -180..-1 deg  -> physical (|theta|, phi+180)
%       main         :    0..180 deg -> physical (theta, phi)
%       bottom mirror: 181..360 deg  -> physical (360-theta, phi+180)
%   Phi extended (3*N_phi cols): [phi-360, phi, phi+360].
%
%   Outputs:
%       theta_ext_deg  : (3N_theta-2 x 1) extended elevation angles [deg].
%       phi_ext_deg    : (3N_phi    x 1) extended azimuth angles [deg].
%       theta_phys_idx : (3N_theta-2 x 1) 1-BASED physical theta index per ext row.
%       phi_offset     : (3N_theta-2 x 1) phi column offset (0 or N_phi/2) per ext row.
%       sin_theta_ext  : (3N_theta-2 x 1) |sin(theta_ext)| solid-angle weight.
%
%   Ported from src/cost/cost_function.py (extended-grid block in
%   build_cost_function / build_directive_physical_masks).
%   Part of: Antenna Array Pattern Optimization Tool.

theta_deg = theta_deg(:);
phi_deg   = phi_deg(:);
N_theta   = numel(theta_deg);
N_phi     = numel(phi_deg);
PHI_180_SHIFT = floor(N_phi / 2);   % [PORT] N_phi // 2

% Extended angle vectors.
% [PORT] -theta_deg[1:][::-1] -> -flip(theta_deg(2:end));
%        360 - theta_deg[N-2::-1] -> 360 - flip(theta_deg(1:end-1)).
theta_ext_deg = [ -flip(theta_deg(2:end)); ...
                   theta_deg; ...
                   360 - flip(theta_deg(1:end-1)) ];
phi_ext_deg   = [ phi_deg - 360; phi_deg; phi_deg + 360 ];

% Physical theta row index for each extended row (1-based MATLAB indices).
% [PORT] [N_theta:-1:2, 1:N_theta, N_theta-1:-1:1]  (already 1-based).
theta_phys_idx = [ (N_theta:-1:2).'; (1:N_theta).'; (N_theta-1:-1:1).' ];

% Phi column offset: mirror rows shift phi by 180deg to reach the back hemisphere.
phi_offset = [ repmat(PHI_180_SHIFT, N_theta - 1, 1); ...
               zeros(N_theta, 1); ...
               repmat(PHI_180_SHIFT, N_theta - 1, 1) ];

% |sin(theta)| on the extended grid (abs keeps mirror rows positive).
sin_theta_ext = abs(sin(deg2rad(theta_ext_deg)));
end
