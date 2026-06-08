function element_patterns = load_element_patterns(directory_path)
% LOAD_ELEMENT_PATTERNS  Load and parse all CST element files from a directory.
%
%   element_patterns = LOAD_ELEMENT_PATTERNS(directory_path)
%
%   Discovers every .txt file in the directory, sorts them by the element index
%   encoded in the filename (e.g. [1], [2], ...), and parses each with
%   parse_cst_file. All elements must share the same (N_theta x N_phi) grid.
%
%   Inputs:
%       directory_path : path to the folder of CST export .txt files.
%
%   Outputs:
%       element_patterns : struct array (1 x N_elements) of parse_cst_file
%                          results, ascending by element index.
%
%   Ported from src/io/cst_parser.py:load_element_patterns.
%   Part of: Antenna Array Pattern Optimization Tool.

listing = dir(fullfile(directory_path, '*.txt'));
% [PORT] Python sorted(dir.glob('*.txt')) + dir() -> MATLAB dir() listing.

if isempty(listing)
    error('load_element_patterns:NoFiles', ...
        'No .txt files found in directory: %s', directory_path);
end

% Sort by the element index extracted from each filename.
n_files   = numel(listing);
elem_idx  = zeros(n_files, 1);
for i = 1:n_files
    elem_idx(i) = extract_element_index(listing(i).name);
end
[~, order] = sort(elem_idx);
listing    = listing(order);

reference_shape = [];
element_patterns = struct([]);
for i = 1:n_files
    txt_file = fullfile(listing(i).folder, listing(i).name);
    pattern  = parse_cst_file(txt_file);

    grid_shape = size(pattern.E_complex);
    if isempty(reference_shape)
        reference_shape = grid_shape;
    elseif ~isequal(grid_shape, reference_shape)
        error('load_element_patterns:ShapeMismatch', ...
            ['Grid shape mismatch across elements: ''%s'' has shape [%s], ' ...
             'expected [%s]. All elements must share the same grid.'], ...
            listing(i).name, num2str(grid_shape), num2str(reference_shape));
    end

    if isempty(element_patterns)
        element_patterns = pattern;
    else
        element_patterns(i) = pattern;
    end
end
end


function idx = extract_element_index(filename)
% Extract the [N] element index from a CST filename.
% Ported from src/io/cst_parser.py:_extract_element_index.
tok = regexp(filename, '\[(\d+)\]', 'tokens', 'once');
if isempty(tok)
    error('load_element_patterns:NoIndex', ...
        ['Could not extract element index from filename ''%s''. ' ...
         'Expected a pattern like ''[N]'' (e.g. ''[3]'').'], filename);
end
idx = str2double(tok{1});
end
