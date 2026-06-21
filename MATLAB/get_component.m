function comp = get_component(pattern, polarization_name)
% GET_COMPONENT  Look up a polarization component's complex grid by name.
%
%   comp = GET_COMPONENT(pattern, polarization_name)
%
%   Inputs:
%       pattern           : element pattern struct returned by parse_cst_file.
%       polarization_name : component name to look up, e.g. 'Copol', 'Cross',
%                            'Theta', 'Phi'. Matched case-insensitively against
%                            fieldnames(pattern.components).
%
%   Outputs:
%       comp : struct with fields abs / phase / complex for the matched component.
%
%   Ported from src/io/cst_parser.py:get_component.
%   Part of: Antenna Array Pattern Optimization Tool.

names = fieldnames(pattern.components);
for i = 1:numel(names)
    if strcmpi(names{i}, polarization_name)
        comp = pattern.components.(names{i});
        return;
    end
end
error('get_component:NotFound', ...
    'Polarization component ''%s'' not found. Available components: %s.', ...
    polarization_name, strjoin(names, ', '));
end
