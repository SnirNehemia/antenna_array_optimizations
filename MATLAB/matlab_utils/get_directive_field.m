function val = get_directive_field(directive, name, default)
% GET_DIRECTIVE_FIELD  Read an optional directive field with a fallback default.
%
%   val = GET_DIRECTIVE_FIELD(directive, name, default)
%
%   Mirrors Python's ``dict.get(key, default)``: returns the field value when the
%   field is present on the struct, otherwise returns ``default``. Directives are
%   represented as MATLAB structs whose absent optional keys are simply missing
%   fields (so ``isfield`` is the analogue of "key present").
%
%   Inputs:
%       directive : scalar struct describing one beam-shaping directive.
%       name      : char/string field name to look up.
%       default   : value returned when the field is absent.
%
%   Outputs:
%       val : directive.(name) when present, else default.
%
%   Part of: Antenna Array Pattern Optimization Tool (shared utility).

% [PORT] Python directive.get(name, default).
if isfield(directive, name)
    val = directive.(name);
else
    val = default;
end
end
