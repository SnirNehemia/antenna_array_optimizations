function idx = nearest_index(grid_deg, target_deg)
% NEAREST_INDEX  Index of the grid point closest to a target angle.
%
%   idx = NEAREST_INDEX(grid_deg, target_deg)
%
%   Inputs:
%       grid_deg   : 1-D angle grid. Units: degrees.
%       target_deg : target angle. Units: degrees.
%
%   Outputs:
%       idx : 1-based index of the nearest grid point (MATLAB convention).
%             Note: the Python original returns the 0-based argmin, so a Python
%             index p corresponds to MATLAB idx = p + 1.
%
%   Ported from src/metrics/metrics.py:_nearest_index (and plotter.py:_nearest_index).
%   Part of: Antenna Array Pattern Optimization Tool.

% [PORT] Python int(np.argmin(abs(grid - target))) -> MATLAB [~, idx] = min(...).
[~, idx] = min(abs(grid_deg(:) - target_deg));
end
