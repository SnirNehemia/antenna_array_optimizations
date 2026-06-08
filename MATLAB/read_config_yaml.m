function config = read_config_yaml(config_path)
% READ_CONFIG_YAML  Minimal YAML reader for this project's config files.
%
%   config = READ_CONFIG_YAML(config_path)
%
%   Base-MATLAB substitute for PyYAML (no YAML toolbox installed). Supports the
%   block-structured subset used by config.yaml / test_config.yaml:
%       - key: value scalars with type inference (int, float incl. 1.0e-6,
%         true/false -> logical, null/~/empty -> [], quoted/unquoted strings)
%       - inline flow lists [a, b] (numeric -> row vector, else cell)
%       - nested mappings via indentation (e.g. optimizer:, output:)
%       - block sequences of mappings via "- " (directives:, scenarios:),
%         including sequences nested inside a scenario mapping
%       - "#" line/inline comments (quote-aware) and blank lines
%
%   Mappings -> struct; sequences -> cell array of structs (the list-of-dicts
%   analogue used throughout the port). Missing/null values -> [].
%
%   Inputs:
%       config_path : path to the YAML file.
%
%   Outputs:
%       config : struct mirroring the parsed YAML mapping.
%
%   Replaces PyYAML usage in scripts/*.py (yaml.safe_load).
%   Part of: Antenna Array Pattern Optimization Tool.

if ~isfile(config_path)
    error('read_config_yaml:FileNotFound', 'Config file not found: %s', config_path);
end

raw = fileread(config_path);
rawlines = regexp(raw, '\r\n|\r|\n', 'split');

% Preprocess: strip comments + blank lines; record indentation.
content = {};
indent  = [];
for i = 1:numel(rawlines)
    line = strip_comment(rawlines{i});
    if isempty(strtrim(line))
        continue;
    end
    lead = regexp(line, '^ *', 'match', 'once');
    content{end + 1, 1} = strtrim(line); %#ok<AGROW>
    indent(end + 1, 1)  = numel(lead);   %#ok<AGROW>
end

n   = numel(content);
pos = 1;

if n == 0
    config = struct();
    return;
end
config = parse_node();


    % ── Recursive parse over the shared cursor `pos` ───────────────
    function val = parse_node()
        if is_dash(content{pos})
            val = parse_sequence(indent(pos));
        else
            val = parse_mapping(indent(pos));
        end
    end

    function m = parse_mapping(base)
        m = struct();
        while pos <= n && indent(pos) == base && ~is_dash(content{pos})
            line = content{pos};
            ci   = find(line == ':', 1, 'first');
            key  = strtrim(line(1:ci - 1));
            rest = strtrim(line(ci + 1:end));
            if isempty(rest)
                % Nested block (mapping or sequence) deeper than `base`.
                pos = pos + 1;
                if pos <= n && indent(pos) > base
                    m.(key) = parse_node();
                else
                    m.(key) = [];   % key with no value
                end
            else
                m.(key) = parse_scalar(rest);
                pos = pos + 1;
            end
        end
    end

    function s = parse_sequence(base)
        s = {};
        while pos <= n && indent(pos) == base && is_dash(content{pos})
            % Rewrite "- key: ..." into a mapping line at base+2, then parse the
            % item (the dash line plus its deeper continuation lines) as a node.
            content{pos} = strtrim(content{pos}(2:end));
            indent(pos)  = base + 2;
            s{end + 1, 1} = parse_node(); %#ok<AGROW>
        end
    end
end


% ────────────────────────── HELPERS ───────────────────────────────

function tf = is_dash(s)
% True for a block-sequence item marker line ("- ..." or bare "-").
tf = strcmp(s, '-') || (numel(s) >= 2 && strcmp(s(1:2), '- '));
end


function out = strip_comment(line)
% Remove a trailing/whole-line # comment, ignoring # inside quotes.
in_single = false;
in_double = false;
out = line;
for i = 1:numel(line)
    c = line(i);
    if c == '''' && ~in_double
        in_single = ~in_single;
    elseif c == '"' && ~in_single
        in_double = ~in_double;
    elseif c == '#' && ~in_single && ~in_double
        out = line(1:i - 1);
        return;
    end
end
end


function v = parse_scalar(s)
% Infer a YAML scalar / inline list value from its text.
s = strtrim(s);

if isempty(s) || strcmp(s, 'null') || strcmp(s, '~')
    v = [];                              % None
    return;
end

% Inline flow list [a, b, ...].
if startsWith(s, '[') && endsWith(s, ']')
    inner = strtrim(s(2:end - 1));
    if isempty(inner)
        v = [];
        return;
    end
    parts = strsplit(inner, ',');
    vals  = cell(1, numel(parts));
    all_numeric = true;
    for i = 1:numel(parts)
        vals{i} = parse_scalar(strtrim(parts{i}));
        if ~(isnumeric(vals{i}) && isscalar(vals{i}))
            all_numeric = false;
        end
    end
    if all_numeric
        v = cell2mat(vals);              % numeric row vector
    else
        v = vals;
    end
    return;
end

% Quoted string -> char (strip quotes).
if (startsWith(s, '"') && endsWith(s, '"')) || (startsWith(s, '''') && endsWith(s, ''''))
    v = s(2:end - 1);
    return;
end

if strcmp(s, 'true'),  v = true;  return; end
if strcmp(s, 'false'), v = false; return; end

num = str2double(s);
if ~isnan(num)
    v = num;                             % int / float / scientific
    return;
end

v = s;                                   % unquoted string
end
