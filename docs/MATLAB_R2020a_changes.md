# Running the MATLAB tool on MATLAB R2020a — copy-paste guide

This tool was written on a newer MATLAB and uses a few functions/properties that
**do not exist in R2020a**. If you get errors like
`Unrecognized property 'GridLineWidth'`, `Unrecognized property 'WordWrap'`,
`Unrecognized function 'clim'`, or `Name argument must be 'ConvertInfAndNaN'`,
apply the edits below.

- **Minimum MATLAB version after these edits:** R2020a
- **Required toolbox:** Optimization Toolbox (nothing else)
- All files are under the `MATLAB/matlab_utils/` and `MATLAB/scripts/` folders.

Just find each **OLD** line and replace it with the **NEW** line. No thinking required.

---

## Change 1 — `clim` → `caxis`  (4 places)

`clim()` was added in R2022a. Use `caxis()` instead.

### 1a. `MATLAB/matlab_utils/ManualWeightsTuner.m` (colorbar init, ~line 421)
**OLD**
```matlab
            clim(obj.ax, [-ManualWeightsTuner.DYNAMIC_RANGE_DB, 0]);
```
**NEW**
```matlab
            caxis(obj.ax, [-ManualWeightsTuner.DYNAMIC_RANGE_DB, 0]);
```

### 1b. `MATLAB/matlab_utils/ManualWeightsTuner.m` (update_heatmap, ~line 1308)
**OLD**
```matlab
            clim(obj.ax, [clim_min, clim_max]);
```
**NEW**
```matlab
            caxis(obj.ax, [clim_min, clim_max]);
```

### 1c. `MATLAB/matlab_utils/save_all_plots.m` (~line 238)
**OLD**
```matlab
colormap(ax, 'jet'); clim(ax, [-dyn 0]);
```
**NEW**
```matlab
colormap(ax, 'jet'); caxis(ax, [-dyn 0]);
```

### 1d. `MATLAB/matlab_utils/save_pattern_gif.m` (~line 77)
**OLD**
```matlab
    colormap(ax_map, 'jet'); clim(ax_map, [-dyn 0]);
```
**NEW**
```matlab
    colormap(ax_map, 'jet'); caxis(ax_map, [-dyn 0]);
```

---

## Change 2 — `GridLineWidth` → `LineWidth`  (1 place)

The `GridLineWidth` axes property was added in R2022b. Before that, the axes
`LineWidth` property controls grid-line thickness.

### `MATLAB/matlab_utils/ManualWeightsTuner.m` (~line 406)
**OLD**
```matlab
            obj.ax.GridLineWidth = 0.5;
```
**NEW**
```matlab
            obj.ax.LineWidth     = 0.5;
```

---

## Change 3 — `uilabel` `WordWrap` → word-wrap helper  (1 place)

The `WordWrap` property of `uilabel` was added in R2020b. Replace the helper
function so it wraps the text itself.

### `MATLAB/matlab_utils/ManualWeightsTuner.m` (the `mk_desc_pos` helper, near the bottom of the file, ~line 1562)
**OLD**
```matlab
% ── File-level helper: wrapped description label at an absolute position ──
function lbl = mk_desc_pos(parent, x, y, w, h, text)
    lbl = uilabel(parent, 'Text', text, 'FontSize', 9, ...
        'FontColor', [0.5 0.5 0.5], 'WordWrap', 'on', ...
        'VerticalAlignment', 'top', 'Position', [x y w h]);
end
```
**NEW** (replace the whole function above with everything below)
```matlab
% ── File-level helper: wrapped description label at an absolute position ──
function lbl = mk_desc_pos(parent, x, y, w, h, text)
    % uilabel 'WordWrap' is R2020b+. For R2020a we pre-wrap the text into a
    % cellstr (rendered as multiple lines) using an estimated character count.
    font_size = 9;
    wrapped   = local_word_wrap(text, w, font_size);
    lbl = uilabel(parent, 'Text', wrapped, 'FontSize', font_size, ...
        'FontColor', [0.5 0.5 0.5], ...
        'VerticalAlignment', 'top', 'Position', [x y w h]);
end

% ── File-level helper: greedy word-wrap to fit a pixel width ──────────
function lines = local_word_wrap(text, width_px, font_size)
    chars_per_line = max(1, floor(width_px / (0.55 * font_size)));
    words = strsplit(char(text));
    lines = {};
    current = '';
    for wi = 1:numel(words)
        word = words{wi};
        if isempty(current)
            candidate = word;
        else
            candidate = [current ' ' word];
        end
        if numel(candidate) > chars_per_line && ~isempty(current)
            lines{end+1} = current; %#ok<AGROW>
            current = word;
        else
            current = candidate;
        end
    end
    if ~isempty(current)
        lines{end+1} = current; %#ok<AGROW>
    end
    if isempty(lines)
        lines = {''};
    end
end
```

---

## Change 4 — `jsonencode(..., 'PrettyPrint', true)` → pretty-print helper  (1 place)

This is the one that causes `Name argument must be 'ConvertInfAndNaN'.`
The `'PrettyPrint'` option was added in R2021a.

### `MATLAB/matlab_utils/run_optimization.m` (in `save_metrics_json`, ~line 205)
**OLD**
```matlab
fprintf(fid, '%s', jsonencode(m, 'PrettyPrint', true));
```
**NEW**
```matlab
fprintf(fid, '%s', pretty_json(jsonencode(m)));
```

Then **add this helper function** to the same file, `run_optimization.m` (paste it
at the end of the file, after the last `end` of the last function):
```matlab
function s = pretty_json(json_str)
% Minimal JSON pretty-printer: indents a compact jsonencode string with
% two-space nesting. Punctuation inside quoted strings is left untouched.
% Replaces the R2021a+ jsonencode 'PrettyPrint' option for R2020a.
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
```

---

## Change 5 (optional but recommended) — friendly version check

Gives a clear message on older MATLAB instead of a cryptic property error.

### `MATLAB/scripts/optimizer_app.m` (right after the header comment, before the path setup)
**ADD** these lines:
```matlab
if verLessThan('matlab', '9.8')
    error('optimizer_app:UnsupportedMATLAB', ...
        'This app requires MATLAB R2020a (9.8) or newer. Detected: %s.', version);
end
```

---

## After applying the changes

1. Start MATLAB R2020a.
2. Make sure the **Optimization Toolbox** is installed (`ver` lists it).
3. `cd` to the `MATLAB/scripts` folder and run:
   ```matlab
   optimizer_app
   ```

That's it — line counts (~421, ~1308, etc.) are approximate; search for the OLD
text if the line number doesn't match exactly.
