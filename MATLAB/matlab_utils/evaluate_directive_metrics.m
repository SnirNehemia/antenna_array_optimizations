function out = evaluate_directive_metrics(af_norm_mag_vm, theta_deg, phi_deg, directive)
% EVALUATE_DIRECTIVE_METRICS  min / SA-weighted-mean / max level in one directive window.
%
%   Uses angular_window_mask on the PHYSICAL grid (no pole/wrap extension) — the
%   compare_classical convention, distinct from evaluate_metrics. The peak_dbi
%   field is the global |AF| peak in dBi; min/mean/max are relative to it.
%
%   Inputs:
%       af_norm_mag_vm : |AF| from power-normalised weights, (N_theta x N_phi) [V/m].
%       theta_deg, phi_deg : angle grids. Units: degrees.
%       directive      : struct with type, theta, optional phi/width/theta_width/phi_width.
%
%   Outputs:
%       out : struct with directive_type, peak_dbi, min_db, mean_db, max_db.
%             min_db/mean_db/max_db are NaN if the window covers no grid points.
%
%   Ported from scripts/compare_classical.py:evaluate_directive_metrics.
%   Part of: Antenna Array Pattern Optimization Tool.

directive_type = directive.type;
target_theta   = directive.theta;
target_phi     = get_directive_field(directive, 'phi', 0.0);
sym_width      = get_directive_field(directive, 'width', []);
theta_width    = get_directive_field(directive, 'theta_width', sym_width);
phi_width      = get_directive_field(directive, 'phi_width',   sym_width);

mask = angular_window_mask(theta_deg, phi_deg, target_theta, target_phi, theta_width, phi_width);

global_peak_vm = max(af_norm_mag_vm, [], 'all');
peak_dbi       = 20 * log10(max(global_peak_vm, 1e-30));

window_values_vm = af_norm_mag_vm(mask);

if isempty(window_values_vm)
    out = struct('directive_type', directive_type, 'peak_dbi', peak_dbi, ...
                 'min_db', NaN, 'mean_db', NaN, 'max_db', NaN);
    return;
end

min_vm = min(window_values_vm);
max_vm = max(window_values_vm);
min_db = 20 * log10(max(min_vm, 1e-30)) - peak_dbi;
max_db = 20 * log10(max(max_vm, 1e-30)) - peak_dbi;

% Solid-angle-weighted mean |AF| in the window.
sin_theta       = sin(deg2rad(theta_deg(:)));         % (N_theta x 1)
weights_sa      = double(mask) .* sin_theta;          % (N_theta x N_phi)
total_weight_sa = sum(weights_sa, 'all');
mean_vm         = sum(af_norm_mag_vm .* weights_sa, 'all') / max(total_weight_sa, 1e-30);
mean_db         = 20 * log10(max(mean_vm, 1e-30)) - peak_dbi;

out = struct('directive_type', directive_type, 'peak_dbi', peak_dbi, ...
             'min_db', min_db, 'mean_db', mean_db, 'max_db', max_db);
end
