function dir_s_db = compute_directivity_trace(run_log, stack1, stack2, theta_deg, phi_deg)
% COMPUTE_DIRECTIVITY_TRACE  Directivity [dBi] toward theta_s at every step of a run.
%
%   DIR_S_DB = COMPUTE_DIRECTIVITY_TRACE(run_log, stack1, stack2, theta_deg, phi_deg)
%
%   Uses the same spherical-integral normalizer as the pattern heatmap
%   (compute_directivity_dbi_grid), so the trace value equals the heatmap
%   value at the steered direction. Shared by save_run_gif and
%   save_run_trace_png so both draw an identical trace.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

n_theta = numel(theta_deg);
n_phi   = numel(phi_deg);
flat1   = reshape(stack1, size(stack1, 1), []);
flat2   = [];
if ~isempty(stack2), flat2 = reshape(stack2, size(stack2, 1), []); end

% Radiated-power Gram G: P_total(w) = w' G w = integral |AF|^2 sin(theta) dOmega.
e_s = run_log.grid.e_s;
sin_theta = sin(deg2rad(theta_deg(:)));
w_sa = repmat(sin_theta, n_phi, 1).';            % 1 x (n_theta*n_phi), theta-fastest
if n_theta > 1, dth = deg2rad(mean(diff(theta_deg))); else, dth = pi;    end
if n_phi   > 1, dph = deg2rad(mean(diff(phi_deg)));   else, dph = 2 * pi; end
G = (flat1 .* w_sa) * flat1';
if ~isempty(flat2), G = G + (flat2 .* w_sa) * flat2'; end
G = G * dth * dph;

T = size(run_log.W, 2);
dir_s_db = zeros(1, T);
for k = 1:T
    w = run_log.W(:, k);
    dir_s_db(k) = 10 * log10(4 * pi * sum(abs(w' * e_s).^2) / real(w' * G * w));
end
end
