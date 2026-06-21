function mask = angular_window_mask(theta_deg, phi_deg, ...
                                    target_theta_deg, target_phi_deg, ...
                                    theta_width_deg, phi_width_deg)
% ANGULAR_WINDOW_MASK  Boolean mask of a rectangular (theta, phi) angular window.
%
%   mask = ANGULAR_WINDOW_MASK(theta_deg, phi_deg, target_theta_deg, ...
%                              target_phi_deg, theta_width_deg, phi_width_deg)
%
%   The window is centred on (target_theta, target_phi) with independent
%   elevation and azimuth half-widths. Operates on whatever theta/phi vectors are
%   passed in; callers needing pole-crossing / 0deg-360deg wrap-around should pass
%   the extended grids built by build_cost_function / build_directive_physical_masks.
%
%   Inputs:
%       theta_deg        : elevation grid, length N_theta. Units: degrees.
%       phi_deg          : azimuth grid,   length N_phi.   Units: degrees.
%       target_theta_deg : window centre elevation. Units: degrees.
%       target_phi_deg   : window centre azimuth.   Units: degrees.
%       theta_width_deg  : full elevation window width. Units: degrees.
%       phi_width_deg    : full azimuth window width.   Units: degrees.
%
%   Outputs:
%       mask : logical (N_theta x N_phi). True where the grid point is in window.
%
%   Ported from src/cost/cost_function.py:angular_window_mask.
%   Part of: Antenna Array Pattern Optimization Tool.

theta_hw = theta_width_deg / 2;
phi_hw   = phi_width_deg   / 2;

% Which theta and phi values fall within the window (column vectors).
theta_in_window = abs(theta_deg(:) - target_theta_deg) <= theta_hw;   % (N_theta x 1)
phi_in_window   = abs(phi_deg(:)   - target_phi_deg)   <= phi_hw;     % (N_phi x 1)

% Outer product: mask(i, j) = theta_i in window AND phi_j in window.
% [PORT] Python theta_in_window[:,None] & phi_in_window[None,:].
mask = theta_in_window & phi_in_window.';   % (N_theta x N_phi)
end
