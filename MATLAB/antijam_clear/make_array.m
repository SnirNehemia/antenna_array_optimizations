function array = make_array(array_folder, component_name, signal_theta_deg, signal_phi_deg)
% ══════════════════════════════════════════════════════════════════
% MAKE_ARRAY
% The antenna and its target, as one object.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   array = MAKE_ARRAY(array_folder, component_name, signal_theta_deg, ...
%                      signal_phi_deg)
%
%   Loads the element patterns and measures the array against the direction it
%   is pointed at, returning both as a single struct. The patterns, the two
%   angle grids and the profile are needed together by everything downstream, so
%   they travel together: every algorithm call in the main script takes 'array'
%   and two or three other things, rather than six.
%
%   THE ONE EXCEPTION, stated so it is a rule and not an inconsistency:
%   steering_vector and array_profile take the patterns and grids explicitly
%   rather than this struct, because this function calls them WHILE building it.
%   Those two are the "measure the antenna" primitives; everything above them
%   takes the bundle.
%
%   Note what is NOT in here: nothing about the jammer, and nothing about the
%   scenario. This struct is knowledge the receiver legitimately has -- its own
%   geometry and where it is trying to listen. It can safely be handed to any
%   algorithm.
%
%   Inputs:
%       array_folder     : folder of CST export .txt files, one per element.
%       component_name   : polarization component, e.g. 'Copol' or 'Theta'.
%                          Required -- see load_array.
%       signal_theta_deg : wanted-signal elevation. Units: degrees.
%       signal_phi_deg   : wanted-signal azimuth.   Units: degrees.
%
%   Outputs:
%       array : struct with fields
%           name               : short array name (the folder's own name), used
%                                in printouts, titles and result filenames.
%           folder             : the full source path, for reporting.
%           component          : the polarization component loaded.
%           element_patterns   : (n_elements x n_theta x n_phi) complex. V/m.
%           theta_deg          : (n_theta x 1) elevation grid. Units: degrees.
%           phi_deg            : (n_phi x 1)   azimuth grid.   Units: degrees.
%           n_elements         : element count.
%           snapshots_per_step : K, snapshots drawn per adaptation step.
%           profile            : struct from array_profile.

% ────────────────────────── CONSTANTS ─────────────────────────────

% Snapshots per element. A covariance estimated from K samples of an n-element
% array lands within about 3 dB of the ideal once K is roughly 2n -- the
% standard Reed-Mallett-Brennan result. K is therefore derived from the array,
% not chosen: 32 snapshots for a 16-element array.
SNAPSHOTS_PER_ELEMENT = 2;

% ────────────────────────── LOAD AND MEASURE ──────────────────────

[element_patterns, theta_deg, phi_deg] = load_array(array_folder, component_name);

n_elements = size(element_patterns, 1);

% Short name for printouts and filenames. The full path is kept too, but a
% figure title or a CSV filename built from an absolute path is unreadable --
% and a Windows path also breaks MATLAB's TeX interpreter on titles.
folder_parts = strsplit(strrep(array_folder, '\', '/'), '/');
folder_parts = folder_parts(~cellfun(@isempty, folder_parts));
array_name   = folder_parts{end};

% Measured before the struct is complete -- hence the explicit arguments.
profile = array_profile(element_patterns, theta_deg, phi_deg, ...
                        signal_theta_deg, signal_phi_deg);

array = struct( ...
    'name',               array_name, ...
    'folder',             array_folder, ...
    'component',          component_name, ...
    'element_patterns',   element_patterns, ...
    'theta_deg',          theta_deg, ...
    'phi_deg',            phi_deg, ...
    'n_elements',         n_elements, ...
    'snapshots_per_step', SNAPSHOTS_PER_ELEMENT * n_elements, ...
    'profile',            profile);
end
