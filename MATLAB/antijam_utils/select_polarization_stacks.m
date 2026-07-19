function [stack1, stack2, label] = select_polarization_stacks(patterns, config)
% SELECT_POLARIZATION_STACKS  Element-pattern stacks per the polarization config.
%
%   [stack1, stack2, label] = SELECT_POLARIZATION_STACKS(patterns, config)
%
%   Mirrors run_optimization's polarization convention: a named component
%   (matched case-insensitively against the components detected in the CST
%   export) or 'total' (incoherent sum of exactly two detected components).
%   Shared by run_antijam and run_jammer_demo.
%
%   Inputs:
%       patterns : struct array from load_element_patterns.
%       config   : parsed config. Required: polarization.
%
%   Outputs:
%       stack1 : (N_el, N_theta, N_phi) complex primary-component stack.
%       stack2 : secondary stack for 'total', else [].
%       label  : printable polarization description.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

if ~isfield(config, 'polarization') || isempty(config.polarization)
    error('select_polarization_stacks:MissingKey', ...
        'Missing required config key: ''polarization''.');
end
pol   = config.polarization;
names = sort(fieldnames(patterns(1).components));
if strcmpi(pol, 'total')
    if numel(names) ~= 2
        error('select_polarization_stacks:BadPolarization', ...
            'polarization ''total'' requires exactly 2 components, found {%s}.', ...
            strjoin(names, ', '));
    end
    stack1 = stack_component(patterns, names{1});
    stack2 = stack_component(patterns, names{2});
    label  = sprintf('total (%s + %s)', names{1}, names{2});
else
    hit = names(strcmpi(names, pol));
    if isempty(hit)
        error('select_polarization_stacks:BadPolarization', ...
            'polarization ''%s'' not found among components {%s}.', ...
            pol, strjoin(names, ', '));
    end
    stack1 = stack_component(patterns, hit{1});
    stack2 = [];
    label  = hit{1};
end
end
