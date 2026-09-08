function prof = kpi_array_profile(stack1, stack2, theta_deg, phi_deg, ...
                                  theta_s_deg, phi_s_deg, loading_db)
% KPI_ARRAY_PROFILE  What a given array can do toward a given target direction.
%
%   prof = KPI_ARRAY_PROFILE(stack1, stack2, theta_deg, phi_deg, ...
%                            theta_s_deg, phi_s_deg, loading_db)
%
%   Measures the properties that the anti-jam stack currently takes as
%   hand-tuned constants, so they can be DERIVED per (array, target) instead.
%   The milestone's constants were each calibrated on one array at one target
%   and do not generalise:
%
%     * `guard_deg` = 5 was derived from the ManyDipoles cut's 15 deg HPBW.
%       Measured across the seven arrays in data/, the beamwidth toward a
%       target spans 24-120 deg, so a single global guard is meaningless -- at
%       the wide end it declares a jammer "outside the main beam" while it sits
%       squarely inside it.
%     * MUSIC's model order is hardcoded `n_sig = 2*n_comp`, which is
%       infeasible whenever the aperture has n_el <= 2*n_comp. `feasible_music`
%       reports that up front so a caller can choose the reactive path rather
%       than discover it mid-run.
%     * The theta-mirror degeneracy (see adapt_cv_init) is a property of the
%       geometry and decides whether a DoA can be interpreted at all.
%
%   Inputs:
%       stack1      : (N_el x N_theta x N_phi) primary far-field stack.
%       stack2      : (N_el x N_theta x N_phi) secondary component, or [].
%       theta_deg   : (1 x N_theta) elevation grid [deg].
%       phi_deg     : (1 x N_phi) azimuth grid [deg].
%       theta_s_deg : target elevation [deg].
%       phi_s_deg   : target azimuth [deg].
%       loading_db  : diagonal loading used to form the quiescent beam [dB].
%
%   Outputs:
%       prof : struct with fields
%           n_el, n_comp        : aperture size and polarization components.
%           dir_s_dbi           : quiescent-beam directivity toward the target
%                                 [dBi] -- what this array can offer HERE.
%           hpbw_theta_deg      : -3 dB beamwidth through the target in theta.
%           hpbw_phi_deg        : the same in phi (NaN if the grid is a cut).
%           guard_deg           : DERIVED exclusion sector = half the wider
%                                 HPBW, floored at one grid step. A jammer
%                                 closer than this is inside the main beam and
%                                 is out of scope by the milestone's own
%                                 definition, whatever the config says.
%           mirror_coh          : coherence with the theta-mirror direction.
%                                 NaN at theta = 90, which is its own mirror.
%           feasible_music      : logical, n_el > 2*n_comp.
%           grid_step_theta_deg, grid_step_phi_deg.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [O].

n_el   = size(stack1, 1);
n_th   = numel(theta_deg);
n_ph   = numel(phi_deg);
n_comp = 1 + ~isempty(stack2);

prof = struct();
prof.n_el   = n_el;
prof.n_comp = n_comp;
prof.feasible_music = n_el > 2 * n_comp;
prof.grid_step_theta_deg = mean(diff(theta_deg));
if n_ph > 1
    prof.grid_step_phi_deg = mean(diff(phi_deg));
else
    prof.grid_step_phi_deg = NaN;
end

% ────────────────────────── QUIESCENT BEAM ────────────────────────
[it_s, ip_s] = nearest_index_2d(theta_deg, phi_deg, theta_s_deg, phi_s_deg);
E1 = reshape(stack1, n_el, []);
E2 = [];
if ~isempty(stack2), E2 = reshape(stack2, n_el, []); end
idx_s = (ip_s - 1) * n_th + it_s;
e_s = E1(:, idx_s);
if ~isempty(E2), e_s = [e_s, E2(:, idx_s)]; end

w = adapt_lcmv_null(eye(n_el), e_s, [], 10^(loading_db / 10));
P = abs(w' * E1).^2;
if ~isempty(E2), P = P + abs(w' * E2).^2; end
P = reshape(P, n_th, n_ph);

% Solid-angle-weighted mean is the isotropic reference, per project convention.
sw = sind(theta_deg(:)) * ones(1, n_ph);
P_avg = sum(P(:) .* sw(:)) / sum(sw(:));
prof.dir_s_dbi = 10 * log10(max(P(it_s, ip_s), realmin) / max(P_avg, realmin));

% ────────────────────────── BEAMWIDTH ─────────────────────────────
prof.hpbw_theta_deg = half_power_width(P(:, ip_s), theta_deg, it_s);
if n_ph > 1
    prof.hpbw_phi_deg = half_power_width(P(it_s, :), phi_deg, ip_s);
else
    prof.hpbw_phi_deg = NaN;
end

% The guard is HALF the wider beamwidth: that is the angular distance from
% boresight at which a source stops being "in the main beam". Floored at one
% grid step so it is always representable on this array's own grid.
wider = max([prof.hpbw_theta_deg, prof.hpbw_phi_deg]);
if ~isfinite(wider)
    wider = prof.hpbw_theta_deg;
end
prof.guard_deg = max(wider / 2, ...
    max(prof.grid_step_theta_deg, 1));

% ────────────────────────── MIRROR DEGENERACY ─────────────────────
% theta = 90 is its own mirror, so the comparison is vacuous there.
if abs(theta_s_deg - 90) < 0.5 * max(prof.grid_step_theta_deg, 1)
    prof.mirror_coh = NaN;
else
    prof.mirror_coh = kpi_steering_coherence(stack1, stack2, theta_deg, phi_deg, ...
        theta_s_deg, phi_s_deg, 180 - theta_s_deg, phi_s_deg);
end
end


% ────────────────────────── HELPERS ───────────────────────────────

function width = half_power_width(cut, axis_deg, i_peak)
% -3 dB width of a 1-D power cut around the sample at i_peak. Returns the full
% span of the axis when the pattern never falls below half power (a very broad
% beam, which is itself the finding on small apertures).
cut = cut(:);
axis_deg = axis_deg(:);
n = numel(cut);
if n < 3 || ~isfinite(cut(i_peak)) || cut(i_peak) <= 0
    width = NaN;
    return
end
half = cut(i_peak) / 2;

lo = find(cut(1:i_peak) < half, 1, 'last');
if isempty(lo), lo = 1; end
hi = i_peak - 1 + find(cut(i_peak:end) < half, 1, 'first');
if isempty(hi), hi = n; end
width = abs(axis_deg(hi) - axis_deg(lo));
end
