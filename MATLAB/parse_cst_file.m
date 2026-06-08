function pattern = parse_cst_file(filepath)
% PARSE_CST_FILE  Parse a single CST Studio far-field ASCII export .txt file.
%
%   pattern = PARSE_CST_FILE(filepath)
%
%   Reads an 8-column whitespace-separated CST export, auto-detects the angular
%   grid resolution, and returns all field components reshaped into
%   (N_theta x N_phi) grids. The complex element pattern is built from the
%   co-polarization magnitude and phase:
%       E_complex(theta, phi) = Abs(Copol) * exp(1i * Phase(Copol) [rad])
%
%   Inputs:
%       filepath : path to the CST far-field export .txt file.
%
%   Outputs:
%       pattern : struct with fields
%           theta_deg     (N_theta x 1)        elevation grid [deg]
%           phi_deg       (N_phi   x 1)        azimuth grid   [deg]
%           E_abs         (N_theta x N_phi)    total E-field magnitude [V/m]
%           copol_abs     (N_theta x N_phi)    co-pol magnitude [V/m]
%           copol_phase   (N_theta x N_phi)    co-pol phase [deg]
%           E_complex     (N_theta x N_phi)    complex co-pol field [V/m]
%           cross_abs     (N_theta x N_phi)    cross-pol magnitude [V/m]
%           cross_phase   (N_theta x N_phi)    cross-pol phase [deg]
%           cross_complex (N_theta x N_phi)    complex cross-pol field [V/m]
%
%   Ported from src/io/cst_parser.py:parse_cst_file.
%   Part of: Antenna Array Pattern Optimization Tool.

% ── Constants (1-based column indices; Python COL_* + 1) ──────────
N_HEADER_LINES  = 2;   % line 1: column labels, line 2: dash separator
N_DATA_COLUMNS  = 8;
COL_THETA       = 1;
COL_PHI         = 2;
COL_E_ABS       = 3;
COL_CROSS_ABS   = 4;
COL_CROSS_PHASE = 5;
COL_COPOL_ABS   = 6;
COL_COPOL_PHASE = 7;
% COL_AXIAL_RATIO = 8 (parsed but unused)

if ~isfile(filepath)
    error('parse_cst_file:FileNotFound', 'CST export file not found: %s', filepath);
end

% Load all numeric data, skipping the 2-line header.
% [PORT] Python np.loadtxt(skiprows=2) -> readmatrix(..., 'NumHeaderLines', 2).
raw_data = readmatrix(filepath, 'NumHeaderLines', N_HEADER_LINES);

if ndims(raw_data) ~= 2 || size(raw_data, 2) ~= N_DATA_COLUMNS %#ok<ISMAT>
    error('parse_cst_file:BadColumns', ...
        'Expected %d columns in ''%s'', got %d.', ...
        N_DATA_COLUMNS, filepath, size(raw_data, 2));
end

% Extract flat angle and field columns.
theta_deg_flat       = raw_data(:, COL_THETA);
phi_deg_flat         = raw_data(:, COL_PHI);
e_abs_flat           = raw_data(:, COL_E_ABS);
cross_abs_flat       = raw_data(:, COL_CROSS_ABS);
cross_phase_deg_flat = raw_data(:, COL_CROSS_PHASE);
copol_abs_flat       = raw_data(:, COL_COPOL_ABS);
copol_phase_deg_flat = raw_data(:, COL_COPOL_PHASE);

[n_theta, n_phi] = detect_grid_shape(theta_deg_flat, phi_deg_flat);

% Sorted 1-D angle vectors from the unique values.
theta_deg = unique(theta_deg_flat);   % (N_theta x 1) ascending
phi_deg   = unique(phi_deg_flat);     % (N_phi   x 1) ascending

% Reshape flat columns into (N_theta x N_phi) grids.
% Theta is the fast (inner) axis in the file, so column-major reshape with
% N_theta first reproduces grid(theta, phi) directly.
% [PORT] Python flat.reshape(n_phi, n_theta).T == MATLAB reshape(flat, n_theta, n_phi).
e_abs_grid           = reshape(e_abs_flat,           n_theta, n_phi);
copol_abs_grid       = reshape(copol_abs_flat,       n_theta, n_phi);
copol_phase_deg_grid = reshape(copol_phase_deg_flat, n_theta, n_phi);
cross_abs_grid       = reshape(cross_abs_flat,       n_theta, n_phi);
cross_phase_deg_grid = reshape(cross_phase_deg_flat, n_theta, n_phi);

% Complex co-pol pattern: convert phase deg -> rad for internal arithmetic.
copol_phase_rad_grid = copol_phase_deg_grid .* (pi / 180);
e_complex_grid       = copol_abs_grid .* exp(1i .* copol_phase_rad_grid);

% Complex cross-pol pattern, same convention.
cross_phase_rad_grid = cross_phase_deg_grid .* (pi / 180);
cross_complex_grid   = cross_abs_grid .* exp(1i .* cross_phase_rad_grid);

pattern = struct( ...
    'theta_deg',     theta_deg, ...
    'phi_deg',       phi_deg, ...
    'E_abs',         e_abs_grid, ...
    'copol_abs',     copol_abs_grid, ...
    'copol_phase',   copol_phase_deg_grid, ...
    'E_complex',     e_complex_grid, ...
    'cross_abs',     cross_abs_grid, ...
    'cross_phase',   cross_phase_deg_grid, ...
    'cross_complex', cross_complex_grid);
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
