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
% Relative power below which a second component is residue rather than a field.
% The measured arrays split into two groups with a 91 dB gap between them
% (-2.4 .. -15.6 dB for genuinely dual-polarized, -107 / -211 dB for ideal
% dipoles), so any threshold in that gap works; -60 dB sits well inside it.
DEGENERATE_COMPONENT_DB = -60.0;

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
    % [P12] Reject a DEGENERATE second component. An ideal-dipole export carries
    % no E_phi, so CST writes numerical residue there — 'data/ManyDipoles' is
    % -107 dB and 'data/Dipole' -211 dB relative to their real component, while
    % every genuinely dual-polarized array here sits between -2.4 and -15.6 dB.
    % Accepting it looks harmless and is not: n_comp becomes 2 while each source
    % is physically rank-1, and adapt_music_doa's hardcoded n_sig = 2*n_comp
    % then reads jammer presence off an eigenvalue that is pure noise. That
    % silently disabled presence detection on 85% of steps and made the
    % predictive nuller fall back to the quiescent beam. No silent default
    % (CLAUDE.md rule 4) — name the component and the fix.
    p1 = sum(abs(stack1(:)).^2);
    p2 = sum(abs(stack2(:)).^2);
    ratio_db = 10 * log10(min(p1, p2) / max(p1, p2));
    if ratio_db < DEGENERATE_COMPONENT_DB
        if p1 < p2, weak = names{1}; strong = names{2};
        else,       weak = names{2}; strong = names{1};
        end
        error('select_polarization_stacks:DegenerateComponent', ...
            ['polarization ''total'' needs two real field components, but ' ...
             '''%s'' carries %.1f dB less power than ''%s'' — it is numerical ' ...
             'residue, not a field, so this array is single-polarization. ' ...
             'Set polarization to ''%s''.'], ...
            weak, -ratio_db, strong, strong);
    end
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
