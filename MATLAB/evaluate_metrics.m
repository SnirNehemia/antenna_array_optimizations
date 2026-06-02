function metrics = evaluate_metrics(element_patterns_stacked, theta_deg, phi_deg, ...
                                    weights_complex, directives, cost_history, ...
                                    precomputed_array_factor, normalizer_power)
% EVALUATE_METRICS  Post-optimization performance metrics for all directives.
%
%   metrics = EVALUATE_METRICS(element_patterns_stacked, theta_deg, phi_deg, ...
%               weights_complex, directives, cost_history, ...
%               precomputed_array_factor, normalizer_power)
%
%   Inputs:
%       element_patterns_stacked  : (N_elements, N_theta, N_phi) complex [V/m].
%                                   Ignored if precomputed_array_factor is given.
%       theta_deg, phi_deg        : angle grids [deg].
%       weights_complex           : optimized complex weights, length N_elements.
%       directives                : cell array of directive structs.
%       cost_history              : cost-per-iteration vector (only its length is used).
%       precomputed_array_factor  : optional (N_theta x N_phi) AF to use directly
%                                   (e.g. sqrt of the total-power grid). Pass [] to
%                                   compute AF from the stack. Default [].
%       normalizer_power          : optional directivity-denominator power. [] uses
%                                   the AF's own radiated power. Default [].
%
%   Outputs:
%       metrics : struct with fields total_cost, iteration_count, global_peak_dbi,
%                 global_peak_theta_deg, global_peak_phi_deg, hpbw_theta_deg,
%                 hpbw_phi_deg, directive_metrics (struct array), peak_to_null_ratio_db.
%                 directive_metrics(k): type, theta_deg, phi_deg, gain_dbi,
%                 null_depth_db ([] for peaks, mirrors Python None), cost_term.
%
%   Ported from src/metrics/metrics.py:evaluate_metrics.
%   Part of: Antenna Array Pattern Optimization Tool.

DEFAULT_TARGET_PHI_DEG   = 0.0;
DEFAULT_DIRECTIVE_WEIGHT = 1.0;

if nargin < 7, precomputed_array_factor = []; end
if nargin < 8, normalizer_power = []; end

theta_deg = theta_deg(:);
phi_deg   = phi_deg(:);

% ── Array factor ──────────────────────────────────────────────────
if ~isempty(precomputed_array_factor)
    array_factor = precomputed_array_factor;
else
    array_factor = compute_array_factor(weights_complex, element_patterns_stacked);
end

power_linear = abs(array_factor) .^ 2;

% ── Directivity grid + global peak ────────────────────────────────
directivity_dbi_grid = compute_directivity_dbi_grid(array_factor, theta_deg, phi_deg, normalizer_power);
global_peak_dbi = max(directivity_dbi_grid, [], 'all');

% Global peak location (1-based -> stored as degrees).
[~, lin_peak]       = max(directivity_dbi_grid(:));
[theta_peak_idx, phi_peak_idx] = ind2sub(size(directivity_dbi_grid), lin_peak);
global_peak_theta_deg = theta_deg(theta_peak_idx);
global_peak_phi_deg   = phi_deg(phi_peak_idx);

n_phi = numel(phi_deg);

% θ-HPBW: cut at peak φ, with pole wrap-around via the opposite azimuth.
if numel(theta_deg) > 1, theta_step = theta_deg(2) - theta_deg(1); else, theta_step = 1.0; end
theta_cut     = directivity_dbi_grid(:, phi_peak_idx);
opp_phi_idx   = mod((phi_peak_idx - 1) + floor(n_phi / 2), n_phi) + 1;   % [PORT] (pp + n_phi//2) % n_phi
opp_theta_cut = directivity_dbi_grid(:, opp_phi_idx);
hpbw_theta_deg = compute_hpbw(theta_cut, theta_peak_idx - 1, opp_theta_cut) * theta_step;

% φ-HPBW: cut at peak θ, rolled so the peak sits at the centre index.
if n_phi > 1, phi_step = phi_deg(2) - phi_deg(1); else, phi_step = 1.0; end
phi_cut_raw = directivity_dbi_grid(theta_peak_idx, :);
roll_shift  = floor(n_phi / 2) - (phi_peak_idx - 1);
phi_cut     = circshift(phi_cut_raw, roll_shift, 2);
hpbw_phi_deg = compute_hpbw(phi_cut, floor(n_phi / 2)) * phi_step;

% ── Per-directive evaluation ──────────────────────────────────────
sin_theta_grid = sin(deg2rad(theta_deg));   % (N_theta x 1), broadcast over phi
phys_masks     = build_directive_physical_masks(theta_deg, phi_deg, directives);

n_dir = numel(directives);
dm = repmat(struct('type', '', 'theta_deg', 0, 'phi_deg', 0, ...
                   'gain_dbi', 0, 'null_depth_db', [], 'cost_term', 0), 1, max(n_dir, 0));
if n_dir == 0
    dm = struct('type', {}, 'theta_deg', {}, 'phi_deg', {}, ...
                'gain_dbi', {}, 'null_depth_db', {}, 'cost_term', {});
end

total_cost     = 0.0;
first_peak_dbi = [];
first_null_dbi = [];

for k = 1:n_dir
    d   = directives{k};
    mask = phys_masks{k};
    directive_type   = d.type;
    target_theta_deg = d.theta;
    target_phi_deg   = get_directive_field(d, 'phi',    DEFAULT_TARGET_PHI_DEG);
    directive_weight = get_directive_field(d, 'weight', DEFAULT_DIRECTIVE_WEIGHT);

    % Realized directivity: strongest point inside the directive window.
    if any(mask(:))
        gain_dbi = max(directivity_dbi_grid(mask));
    else
        theta_idx = nearest_index(theta_deg, target_theta_deg);
        phi_idx   = nearest_index(phi_deg,   target_phi_deg);
        gain_dbi  = directivity_dbi_grid(theta_idx, phi_idx);
    end

    % Solid-angle-weighted mean power over the same window.
    weights_sa   = double(mask) .* sin_theta_grid;
    total_weight = sum(weights_sa, 'all');
    if total_weight > 1e-30
        mean_window_power = sum(power_linear .* weights_sa, 'all') / max(total_weight, 1e-30);
    else
        mean_window_power = 0.0;
    end

    if strcmp(directive_type, 'peak')
        cost_term     = -directive_weight * mean_window_power;
        null_depth_db = [];                       % mirrors Python None
        if isempty(first_peak_dbi), first_peak_dbi = gain_dbi; end
    else   % 'null'
        cost_term     = +directive_weight * mean_window_power;
        null_depth_db = gain_dbi - global_peak_dbi;
        if isempty(first_null_dbi), first_null_dbi = gain_dbi; end
    end

    total_cost = total_cost + cost_term;

    dm(k).type          = directive_type;
    dm(k).theta_deg     = double(target_theta_deg);
    dm(k).phi_deg       = double(target_phi_deg);
    dm(k).gain_dbi      = double(gain_dbi);
    dm(k).null_depth_db = null_depth_db;
    dm(k).cost_term     = double(cost_term);
end

% ── Peak-to-null ratio ────────────────────────────────────────────
if ~isempty(first_peak_dbi) && ~isempty(first_null_dbi)
    peak_to_null_ratio_db = first_peak_dbi - first_null_dbi;
else
    peak_to_null_ratio_db = [];
end

metrics = struct();
metrics.total_cost            = double(total_cost);
metrics.iteration_count       = numel(cost_history);
metrics.global_peak_dbi       = double(global_peak_dbi);
metrics.global_peak_theta_deg = double(global_peak_theta_deg);
metrics.global_peak_phi_deg   = double(global_peak_phi_deg);
metrics.hpbw_theta_deg        = double(hpbw_theta_deg);
metrics.hpbw_phi_deg          = double(hpbw_phi_deg);
metrics.directive_metrics     = dm;
metrics.peak_to_null_ratio_db = peak_to_null_ratio_db;
end
