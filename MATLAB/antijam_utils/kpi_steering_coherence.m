function coh = kpi_steering_coherence(stack1, stack2, theta_deg, phi_deg, ...
                                      theta_s_deg, phi_s_deg, theta_j_deg, phi_j_deg)
% KPI_STEERING_COHERENCE  Normalized coherence between the target and jammer steering columns.
%
%   coh = KPI_STEERING_COHERENCE(stack1, stack2, theta_deg, phi_deg, ...
%                                theta_s_deg, phi_s_deg, theta_j_deg, phi_j_deg)
%
%   [P12] The array-independent difficulty coordinate. Angular separation is not
%   comparable across arrays — 45 deg means something different on a 6-element
%   array than on a 4x4 — but coherence is: it is the quantity that actually
%   sets how hard the jammer is to null without taking the target with it. A
%   coherence SPIKE at an angle where separation says the problem should be easy
%   is the signature of a grating lobe, which is precisely the blind spot the
%   0.6-lambda 4x4 array is suspected of carrying.
%
%       coh = ||e_s' * e_j||_F / (||e_s||_F * ||e_j||_F)
%
%   With one polarization component this is the familiar
%   |e_s' e_j| / (|e_s| |e_j|). With two (polarization 'total', n_comp = 2) the
%   Frobenius form generalizes it to the subspace pair and stays in [0, 1] by
%   Cauchy-Schwarz. 0 = orthogonal (easiest), 1 = indistinguishable (hopeless).
%
%   Steering columns are built exactly as sim_engine_init builds e_s — same
%   flatten order, same nearest_index_2d snap — so this reports the coherence
%   the engine actually sees, not the one the requested angles imply.
%
%   Inputs:
%       stack1, stack2 : (N_el x N_theta x N_phi) far-field stacks; stack2 is []
%                        for single-component operation. Units: V/m.
%       theta_deg, phi_deg : grid axes. Units: degrees.
%       theta_s_deg, phi_s_deg : target direction. Units: degrees.
%       theta_j_deg, phi_j_deg : jammer direction. Units: degrees.
%
%   Outputs:
%       coh : normalized coherence in [0, 1]. Units: dimensionless.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P12].

e_s = steering_column(stack1, stack2, theta_deg, phi_deg, theta_s_deg, phi_s_deg);
e_j = steering_column(stack1, stack2, theta_deg, phi_deg, theta_j_deg, phi_j_deg);

den = norm(e_s, 'fro') * norm(e_j, 'fro');
if den == 0
    error('kpi_steering_coherence:NullPattern', ...
        ['A steering column is identically zero (target %.1f/%.1f, jammer ' ...
         '%.1f/%.1f deg). The element patterns carry no field in that ' ...
         'direction, so coherence is undefined.'], ...
        theta_s_deg, phi_s_deg, theta_j_deg, phi_j_deg);
end
coh = norm(e_s' * e_j, 'fro') / den;
coh = min(coh, 1);          % clamp for floating-point round-off
end


function e = steering_column(stack1, stack2, theta_deg, phi_deg, th_deg, ph_deg)
% One (N_el x n_comp) steering column, mirroring sim_engine_init's e_s build.
n_theta = numel(theta_deg);
n_phi   = numel(phi_deg);
n_el    = size(stack1, 1);
E1 = reshape(stack1, n_el, n_theta * n_phi);

[it, ip] = nearest_index_2d(theta_deg, phi_deg, th_deg, ph_deg);
idx = (ip - 1) * n_theta + it;      % linear index into the (theta,phi) flatten
e   = E1(:, idx);
if ~isempty(stack2)
    E2 = reshape(stack2, n_el, n_theta * n_phi);
    e  = [e, E2(:, idx)];
end
end
