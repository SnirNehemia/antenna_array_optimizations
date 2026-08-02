function save_run_params(scn_cfg, aj, sim_cfg, gif_cfg, methods, pol, out_path)
% SAVE_RUN_PARAMS  Dump one demo run's inputs to a pretty-printed JSON file.
%
%   SAVE_RUN_PARAMS(scn_cfg, aj, sim_cfg, gif_cfg, methods, pol, out_path)
%
%   Records the scenario definition (scn_cfg, as authored in the calling
%   script — before sim_scenario expands it into a per-step timeline), the
%   antijam config (aj, from config.yaml), the closed-loop sim settings
%   (sim_cfg), the video settings (gif_cfg), the algorithm list, and the
%   polarization mode, so a later look at a timestamped results folder shows
%   exactly what produced it without re-reading the script/config as they
%   stood at run time.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

params = struct( ...
    'timestamp',  char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'scenario',   scn_cfg, ...
    'antijam',    aj, ...
    'sim',        sim_cfg, ...
    'video',      gif_cfg, ...
    'methods',    {methods}, ...
    'polarization', pol);

fid = fopen(out_path, 'w');
% [MATLAB] jsonencode's 'PrettyPrint' name-value option is R2021a+. For
% R2020a compatibility, encode compactly and indent with a local helper
% (same pattern as run_optimization.m's save_metrics_json).
fprintf(fid, '%s', pretty_json(jsonencode(params)));
fclose(fid);
end


% ────────────────────────── HELPERS ───────────────────────────────

function s = pretty_json(json_str)
% Minimal JSON pretty-printer: indents a compact jsonencode string with
% two-space nesting. Punctuation inside quoted strings is left untouched.
json_str    = char(json_str);
indent_unit = '  ';
out         = '';
level       = 0;
in_string   = false;
escaped     = false;
n           = numel(json_str);
i           = 1;
while i <= n
    c = json_str(i);
    if in_string
        out(end+1) = c; %#ok<AGROW>
        if escaped
            escaped = false;
        elseif c == '\'
            escaped = true;
        elseif c == '"'
            in_string = false;
        end
        i = i + 1;
        continue;
    end
    switch c
        case '"'
            in_string  = true;
            out(end+1) = c; %#ok<AGROW>
        case {'{', '['}
            % Keep empty containers ({} / []) on one line.
            if i < n && ((c == '{' && json_str(i+1) == '}') || ...
                         (c == '[' && json_str(i+1) == ']'))
                out(end+1) = c;             %#ok<AGROW>
                out(end+1) = json_str(i+1); %#ok<AGROW>
                i = i + 2;
                continue;
            end
            level      = level + 1;
            out(end+1) = c; %#ok<AGROW>
            out = [out newline repmat(indent_unit, 1, level)];
        case {'}', ']'}
            level = level - 1;
            out   = [out newline repmat(indent_unit, 1, level)];
            out(end+1) = c; %#ok<AGROW>
        case ','
            out(end+1) = c; %#ok<AGROW>
            out = [out newline repmat(indent_unit, 1, level)];
        case ':'
            out = [out ': '];
        otherwise
            out(end+1) = c; %#ok<AGROW>
    end
    i = i + 1;
end
s = out;
end
