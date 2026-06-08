function masks = build_directive_physical_masks(theta_deg, phi_deg, directives)
% BUILD_DIRECTIVE_PHYSICAL_MASKS  Physical-grid boolean mask for each directive.
%
%   masks = BUILD_DIRECTIVE_PHYSICAL_MASKS(theta_deg, phi_deg, directives)
%
%   Uses the same extended-grid logic as build_cost_function so masks match the
%   optimizer's angular windows exactly, including phi 0/360 wrap-around and
%   theta pole-crossing at both poles.
%
%   Inputs:
%       theta_deg  : elevation grid, length N_theta. Units: degrees.
%       phi_deg    : azimuth grid,   length N_phi.   Units: degrees.
%       directives : cell array of directive structs (the list-of-dicts analogue).
%
%   Outputs:
%       masks : cell array (numel(directives) x 1) of logical (N_theta x N_phi).
%
%   Ported from src/cost/cost_function.py:build_directive_physical_masks.
%   Part of: Antenna Array Pattern Optimization Tool.

DEFAULT_TARGET_PHI_DEG = 0.0;

theta_deg = theta_deg(:);
phi_deg   = phi_deg(:);
N_theta   = numel(theta_deg);
N_phi     = numel(phi_deg);

[theta_ext_deg, phi_ext_deg, theta_phys_idx, phi_offset] = ...
    extended_grid_maps(theta_deg, phi_deg);

masks = cell(numel(directives), 1);
for k = 1:numel(directives)
    d  = directives{k};
    tt = d.theta;
    tp = get_directive_field(d, 'phi', DEFAULT_TARGET_PHI_DEG);
    sym_width = get_directive_field(d, 'width', []);
    tw = get_directive_field(d, 'theta_width', sym_width);
    pw = get_directive_field(d, 'phi_width',   sym_width);

    ext_mask = angular_window_mask(theta_ext_deg, phi_ext_deg, tt, tp, tw, pw);
    [ext_t, ext_p] = find(ext_mask);     % 1-based row/col indices

    % Map extended indices -> physical grid indices.
    % [PORT] phys_phi = mod(mod(ext_p-1, N_phi) + phi_offset(ext_t), N_phi) + 1.
    phys_theta = theta_phys_idx(ext_t);
    phys_phi   = mod(mod(ext_p - 1, N_phi) + phi_offset(ext_t), N_phi) + 1;

    m = false(N_theta, N_phi);
    m(sub2ind([N_theta, N_phi], phys_theta, phys_phi)) = true;
    masks{k} = m;
end
end
