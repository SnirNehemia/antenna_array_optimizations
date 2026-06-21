function stack = stack_component(patterns, component_name)
% STACK_COMPONENT  Stack one polarization component's complex grid across elements.
%
%   stack = STACK_COMPONENT(patterns, component_name)
%
%   Inputs:
%       patterns       : struct array of element patterns (parse_cst_file output).
%       component_name : component name to look up via get_component, e.g.
%                         'Copol', 'Cross', 'Theta', 'Phi'. Matched
%                         case-insensitively.
%
%   Outputs:
%       stack : (N_elements x N_theta x N_phi) complex array.
%
%   Part of: Antenna Array Pattern Optimization Tool.

n = numel(patterns);
comp1 = get_component(patterns(1), component_name);
[n_theta, n_phi] = size(comp1.complex);
stack = zeros(n, n_theta, n_phi);
stack(1, :, :) = comp1.complex;
for i = 2:n
    comp = get_component(patterns(i), component_name);
    stack(i, :, :) = comp.complex;
end
end
