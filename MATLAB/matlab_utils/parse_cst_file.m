function pattern = parse_cst_file(filepath)
% PARSE_CST_FILE  Parse a single CST Studio far-field ASCII export .txt file.
%
%   pattern = PARSE_CST_FILE(filepath)
%
%   Reads an 8-column whitespace-separated CST export, auto-detects the angular
%   grid resolution, and returns all field components reshaped into
%   (N_theta x N_phi) grids. Polarization components are detected generically
%   from the header's Abs(<name>)/Phase(<name>) column pairs (e.g.
%   Copol/Cross, or Theta/Phi), so this function supports any CST far-field
%   export basis without code changes.
%
%   For each detected component <name>, the complex field is built as:
%       complex(theta, phi) = Abs(<name>) * exp(1i * Phase(<name>) [rad])
%
%   Inputs:
%       filepath : path to the CST far-field export .txt file.
%
%   Outputs:
%       pattern : struct with fields
%           theta_deg   (N_theta x 1)        elevation grid [deg]
%           phi_deg     (N_phi   x 1)        azimuth grid   [deg]
%           E_abs       (N_theta x N_phi)    magnitude-only total field
%                                             (Abs(E) or Abs(Grlz))
%           axial_ratio (N_theta x N_phi)    axial ratio [dimensionless]
%           components  struct, one field per detected Abs(<name>)/
%                       Phase(<name>) pair, keyed by <name> exactly as
%                       written in the header (e.g. Copol, Cross, Theta,
%                       Phi). Each value is a struct with:
%                           abs     (N_theta x N_phi)  magnitude
%                           phase   (N_theta x N_phi)  phase [deg]
%                           complex (N_theta x N_phi)  abs .* exp(1i*phase_rad)
%
%   Ported from src/io/cst_parser.py:parse_cst_file.
%   Part of: Antenna Array Pattern Optimization Tool.

N_HEADER_LINES = 2;   % line 1: column labels, line 2: dash separator
N_DATA_COLUMNS = 8;
COL_THETA      = 1;
COL_PHI        = 2;

if ~isfile(filepath)
    error('parse_cst_file:FileNotFound', 'CST export file not found: %s', filepath);
end

% Read the header line (column labels) to identify column roles generically.
fid = fopen(filepath, 'r');
header_line = fgetl(fid);
fclose(fid);
header_info = parse_header_columns(header_line, COL_PHI);

% Load all numeric data, skipping the 2-line header.
% [PORT] Python np.loadtxt(skiprows=2) -> readmatrix(..., 'NumHeaderLines', 2).
raw_data = readmatrix(filepath, 'NumHeaderLines', N_HEADER_LINES);

if ndims(raw_data) ~= 2 || size(raw_data, 2) ~= N_DATA_COLUMNS %#ok<ISMAT>
    error('parse_cst_file:BadColumns', ...
        'Expected %d columns in ''%s'', got %d.', ...
        N_DATA_COLUMNS, filepath, size(raw_data, 2));
end

theta_deg_flat = raw_data(:, COL_THETA);
phi_deg_flat   = raw_data(:, COL_PHI);

[n_theta, n_phi] = detect_grid_shape(theta_deg_flat, phi_deg_flat);

% Sorted 1-D angle vectors from the unique values.
theta_deg = unique(theta_deg_flat);   % (N_theta x 1) ascending
phi_deg   = unique(phi_deg_flat);     % (N_phi   x 1) ascending

% Reshape flat columns into (N_theta x N_phi) grids.
% Theta is the fast (inner) axis in the file, so column-major reshape with
% N_theta first reproduces grid(theta, phi) directly.
% [PORT] Python flat.reshape(n_phi, n_theta).T == MATLAB reshape(flat, n_theta, n_phi).
e_abs_grid       = reshape(raw_data(:, header_info.e_abs_col),       n_theta, n_phi);
axial_ratio_grid = reshape(raw_data(:, header_info.axial_ratio_col), n_theta, n_phi);

% Build the complex field for every detected polarization component.
% [PORT] phase_rad_grid = phase_deg_grid .* (pi / 180);
% [PORT] complex_grid   = abs_grid .* exp(1i .* phase_rad_grid);
components = struct();
comp_names = fieldnames(header_info.components);
for i = 1:numel(comp_names)
    name = comp_names{i};
    cols = header_info.components.(name);
    abs_grid       = reshape(raw_data(:, cols.abs_col),   n_theta, n_phi);
    phase_deg_grid = reshape(raw_data(:, cols.phase_col), n_theta, n_phi);
    phase_rad_grid = phase_deg_grid .* (pi / 180);
    complex_grid   = abs_grid .* exp(1i .* phase_rad_grid);
    components.(name) = struct( ...
        'abs',     abs_grid, ...
        'phase',   phase_deg_grid, ...
        'complex', complex_grid);
end

pattern = struct( ...
    'theta_deg',   theta_deg, ...
    'phi_deg',     phi_deg, ...
    'E_abs',       e_abs_grid, ...
    'axial_ratio', axial_ratio_grid, ...
    'components',  components);
end


function header_info = parse_header_columns(header_line, col_phi)
% Identify the role of each data column from the CST header line.
%
% Searches the header for Abs(<name>)/Phase(<name>) pairs (e.g.
% Abs(Copol)/Phase(Copol), or Abs(Theta)/Phase(Theta)) and the Ax.Ratio
% column. Columns 1 and 2 are always Theta and Phi. The first Abs(<name>)
% column with no matching Phase(<name>) is treated as the magnitude-only
% total-field column (Abs(E) or Abs(Grlz)).
%
% Returns a struct with fields:
%   e_abs_col       : 1-based column index of the magnitude-only total field.
%   axial_ratio_col : 1-based column index of Ax.Ratio.
%   components      : struct, one field per detected component name, each
%                      holding a struct with abs_col / phase_col indices.
%
% Ported from src/io/cst_parser.py:_parse_header_columns.
field_pattern = 'Abs\(([^)]*)\)|Phase\(([^)]*)\)|Ax\.Ratio';
matches = regexp(header_line, field_pattern, 'match');
tokens  = regexp(header_line, field_pattern, 'tokens');

abs_cols   = struct();   % name -> column index
phase_cols = struct();   % name -> column index
axial_ratio_col = [];

col_idx = col_phi + 1;  % field columns start right after Theta, Phi
for i = 1:numel(matches)
    if strcmp(matches{i}, 'Ax.Ratio')
        axial_ratio_col = col_idx;
    elseif startsWith(matches{i}, 'Abs(')
        name = matlab.lang.makeValidName(strtrim(tokens{i}{1}));
        if isfield(abs_cols, name)
            error('parse_cst_file:DuplicateColumn', ...
                'Duplicate ''Abs(%s)'' column in header: %s', name, header_line);
        end
        abs_cols.(name) = col_idx;
    else
        name = matlab.lang.makeValidName(strtrim(tokens{i}{1}));
        if isfield(phase_cols, name)
            error('parse_cst_file:DuplicateColumn', ...
                'Duplicate ''Phase(%s)'' column in header: %s', name, header_line);
        end
        phase_cols.(name) = col_idx;
    end
    col_idx = col_idx + 1;
end

% The magnitude-only total-field column is the Abs(...) entry with no
% matching Phase(...) entry (e.g. Abs(E), Abs(Grlz)).
abs_names   = fieldnames(abs_cols);
phase_names = fieldnames(phase_cols);
e_abs_names = setdiff(abs_names, phase_names);
if numel(e_abs_names) ~= 1
    error('parse_cst_file:BadHeader', ...
        ['Expected exactly one magnitude-only Abs(...) column (e.g. ''Abs(E)'' ' ...
         'or ''Abs(Grlz)''), found %d in header: %s'], numel(e_abs_names), header_line);
end
e_abs_col = abs_cols.(e_abs_names{1});
abs_cols  = rmfield(abs_cols, e_abs_names{1});

if isempty(axial_ratio_col)
    error('parse_cst_file:BadHeader', 'No ''Ax.Ratio'' column found in header: %s', header_line);
end

components = struct();
remaining_abs_names = fieldnames(abs_cols);
for i = 1:numel(remaining_abs_names)
    name = remaining_abs_names{i};
    if ~isfield(phase_cols, name)
        error('parse_cst_file:BadHeader', ...
            '''Abs(%s)'' column has no matching ''Phase(%s)'' column in header: %s', ...
            name, name, header_line);
    end
    components.(name) = struct('abs_col', abs_cols.(name), 'phase_col', phase_cols.(name));
end

if isempty(fieldnames(components))
    error('parse_cst_file:BadHeader', ...
        'No Abs(<name>)/Phase(<name>) component pairs found in header: %s', header_line);
end

header_info = struct( ...
    'e_abs_col',       e_abs_col, ...
    'axial_ratio_col', axial_ratio_col, ...
    'components',      components);
end


function [n_theta, n_phi] = detect_grid_shape(theta_deg_flat, phi_deg_flat)
% Infer (N_theta, N_phi) from flat angle arrays; error on incomplete grids.
% Ported from src/io/cst_parser.py:_detect_grid_shape.
n_theta = numel(unique(theta_deg_flat));
n_phi   = numel(unique(phi_deg_flat));
n_rows  = numel(theta_deg_flat);
if n_theta * n_phi ~= n_rows
    error('parse_cst_file:GridMismatch', ...
        ['Grid shape mismatch: %d unique Theta x %d unique Phi = %d, ' ...
         'but file has %d data rows. Expected a complete rectangular grid.'], ...
        n_theta, n_phi, n_theta * n_phi, n_rows);
end
end
