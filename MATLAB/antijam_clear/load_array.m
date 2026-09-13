function [element_patterns, theta_deg, phi_deg] = load_array(array_folder, component_name)
% ══════════════════════════════════════════════════════════════════
% LOAD_ARRAY
% Load one polarization component of a CST element-pattern folder.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   [element_patterns, theta_deg, phi_deg] = LOAD_ARRAY(array_folder, component_name)
%
%   Parses every CST .txt export in array_folder (via matlab_utils) and stacks
%   ONE polarization component into a single complex cube. That cube is the only
%   description of the antenna the rest of this folder ever uses: every question
%   any function asks -- "what does element n do in direction (theta, phi)?" --
%   is answered by indexing it.
%
%   These are measured CST patterns, not an analytic exp(1i*k*d*sin(theta))
%   model. The elements are neither identical nor isotropic, and that is the
%   whole point: the results are statements about this hardware.
%
%   POLARIZATION. Only one component is loaded. The arrays here use two export
%   conventions -- the patch arrays carry Copol/Cross, the dipole and monopole
%   arrays carry Theta/Phi -- so there is no correct default and none is
%   supplied (CLAUDE.md rule 4). The caller names the component on a visible
%   line.
%
%   [FUTURE] Two-component ("total") operation is a documented extension, not a
%   hidden option: call this twice and carry a second cube. Each source then has
%   TWO steering vectors instead of one, and every downstream formula gains a
%   sum over components -- see steering_vector for the exact shape change.
%
%   Inputs:
%       array_folder   : path to a folder of CST export .txt files, one per
%                        element (e.g. 'data/Monopoles').
%       component_name : polarization component to load, matched
%                        case-insensitively against the components present in
%                        the export. 'Copol' | 'Cross' | 'Theta' | 'Phi'.
%
%   Outputs:
%       element_patterns : (n_elements x n_theta x n_phi) complex far field.
%                          Units: V/m (arbitrary common scale).
%       theta_deg        : (n_theta x 1) elevation grid. Units: degrees.
%       phi_deg          : (n_phi x 1)   azimuth grid.   Units: degrees.

% ────────────────────────── PARSE ─────────────────────────────────

if nargin < 2 || isempty(component_name)
    error('load_array:MissingComponent', ...
        ['Missing required argument ''component_name''. There is no correct ' ...
         'default: patch arrays export Copol/Cross, dipole and monopole ' ...
         'arrays export Theta/Phi. Name the component explicitly.']);
end
if ~exist(array_folder, 'dir')
    error('load_array:NoSuchFolder', 'Array folder not found: %s', array_folder);
end

patterns = load_element_patterns(array_folder);

% ────────────────────────── STACK ONE COMPONENT ───────────────────

% stack_component raises a descriptive error naming the available components
% if component_name is not one of them, so no check is duplicated here.
element_patterns = stack_component(patterns, component_name);

theta_deg = patterns(1).theta_deg(:);
phi_deg   = patterns(1).phi_deg(:);
end
