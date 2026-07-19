function [E1c, E2c, cut_ang] = extract_cut(stack1, stack2, theta_deg, phi_deg, aj)
% EXTRACT_CUT  The engine's 1-D azimuth cut per antijam.cut_type.
%
%   [E1c, E2c, cut_ang] = EXTRACT_CUT(stack1, stack2, theta_deg, phi_deg, aj)
%
%   'phi_cut'   : azimuth sweep at fixed theta = aj.cut_theta_deg;
%                 cut_ang = phi grid.
%   'theta_cut' : principal plane combining the phi = cut_phi_deg and +180
%                 half-planes into a -90..+90 deg axis (complex analogue of
%                 principal_plane_cut).
%   Shared by run_antijam and run_jammer_demo.
%
%   Inputs:
%       stack1, stack2      : (N_el, N_theta, N_phi) stacks (stack2 may be []).
%       theta_deg, phi_deg  : full angle grids [deg].
%       aj                  : antijam config (cut_type + cut_theta_deg /
%                             cut_phi_deg).
%
%   Outputs:
%       E1c, E2c : (N_el x N_angles) complex cuts (E2c [] when stack2 is []).
%       cut_ang  : (1 x N_angles) cut angle axis [deg].
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

switch aj.cut_type
    case 'phi_cut'
        if ~isfield(aj, 'cut_theta_deg') || isempty(aj.cut_theta_deg)
            error('extract_cut:MissingKey', 'phi_cut requires antijam.cut_theta_deg.');
        end
        it = nearest_index(theta_deg(:), aj.cut_theta_deg);
        cut_ang = phi_deg(:).';
        E1c = squeeze(stack1(:, it, :));
        E2c = [];
        if ~isempty(stack2), E2c = squeeze(stack2(:, it, :)); end
    case 'theta_cut'
        if ~isfield(aj, 'cut_phi_deg') || isempty(aj.cut_phi_deg)
            error('extract_cut:MissingKey', 'theta_cut requires antijam.cut_phi_deg.');
        end
        ip0   = nearest_index(phi_deg(:), mod(aj.cut_phi_deg, 360));
        ip180 = nearest_index(phi_deg(:), mod(aj.cut_phi_deg + 180, 360));
        cut_ang = [-fliplr(theta_deg(2:end)), theta_deg];
        E1c = [fliplr(squeeze(stack1(:, 2:end, ip180))), squeeze(stack1(:, :, ip0))];
        E2c = [];
        if ~isempty(stack2)
            E2c = [fliplr(squeeze(stack2(:, 2:end, ip180))), squeeze(stack2(:, :, ip0))];
        end
    otherwise
        error('extract_cut:BadCutType', 'Unknown cut_type ''%s''.', aj.cut_type);
end
end
