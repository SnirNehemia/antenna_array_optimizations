function cost_fn = build_cost_function(element_patterns_stacked, theta_deg, phi_deg, ...
                                       directives, mode, element_patterns_secondary)
% BUILD_COST_FUNCTION  Build the composite cost J(x) as a function handle.
%
%   cost_fn = BUILD_COST_FUNCTION(element_patterns_stacked, theta_deg, phi_deg, ...
%                                 directives, mode, element_patterns_secondary)
%
%   Constructs the extended (theta, phi) grid (pole mirrors + phi tiling), reduces
%   each directive's angular window to a sparse list of physical grid points, and
%   returns cost_fn(x) -> scalar J. Composite cost:
%       J(x) = sum_k  lambda_k * C_k(x)
%
%   Inputs:
%       element_patterns_stacked   : (N_elements, N_theta, N_phi) complex [V/m].
%       theta_deg                  : (N_theta x 1) elevation grid [deg], 0..180.
%       phi_deg                    : (N_phi   x 1) azimuth grid [deg], 0..<360, N_phi even.
%       directives                 : cell array of directive structs.
%       mode                       : 'standard' | 'phase_only' | 'amplitude_only'.
%       element_patterns_secondary : optional cross-pol stack (same size) for the
%                                    total-power mode |AF_primary|^2 + |AF_secondary|^2.
%                                    Pass [] for single-polarisation. Default [].
%
%   Outputs:
%       cost_fn : function handle cost_fn(x) -> double J.
%
%   Ported from src/cost/cost_function.py:build_cost_function.
%   Part of: Antenna Array Pattern Optimization Tool.

if nargin < 5 || isempty(mode)
    mode = 'standard';
end
if nargin < 6
    element_patterns_secondary = [];
end

VALID_MODES           = {'standard', 'phase_only', 'amplitude_only'};
VALID_DIRECTIVE_TYPES = {'peak', 'null'};
VALID_AGGREGATIONS    = {'mean', 'max', 'min'};
DEFAULT_DIRECTIVE_WEIGHT = 1.0;
DEFAULT_TARGET_PHI_DEG   = 0.0;

if ~ismember(mode, VALID_MODES)
    error('build_cost_function:BadMode', ...
        'Unknown optimization mode ''%s''. Valid: standard, phase_only, amplitude_only.', mode);
end

theta_deg = theta_deg(:);
phi_deg   = phi_deg(:);
N_theta   = numel(theta_deg);
N_phi     = numel(phi_deg);

[theta_ext_deg, phi_ext_deg, theta_phys_idx, phi_offset, sin_theta_ext] = ...
    extended_grid_maps(theta_deg, phi_deg);

% ── Sparse mask pre-computation (one entry per directive) ─────────
mask_data = struct('lin_idx', {}, 'sin_w', {}, 'aggregation', {}, ...
                   'd_type', {}, 'd_weight', {});
for k = 1:numel(directives)
    d  = directives{k};
    tt = d.theta;
    tp = get_directive_field(d, 'phi', DEFAULT_TARGET_PHI_DEG);
    sym_width = get_directive_field(d, 'width', []);
    tw = get_directive_field(d, 'theta_width', sym_width);
    pw = get_directive_field(d, 'phi_width',   sym_width);

    ext_mask = angular_window_mask(theta_ext_deg, phi_ext_deg, tt, tp, tw, pw);
    [ext_t, ext_p] = find(ext_mask);

    phys_theta = theta_phys_idx(ext_t);
    phys_phi   = mod(mod(ext_p - 1, N_phi) + phi_offset(ext_t), N_phi) + 1;
    lin_idx    = sub2ind([N_theta, N_phi], phys_theta, phys_phi);

    aggregation = get_directive_field(d, 'aggregation', 'mean');
    d_type      = d.type;
    d_weight    = get_directive_field(d, 'weight', DEFAULT_DIRECTIVE_WEIGHT);

    % Validate at build time, not inside the hot loop.
    if ~ismember(d_type, VALID_DIRECTIVE_TYPES)
        error('build_cost_function:BadType', ...
            'Unknown directive type ''%s''. Valid: peak, null.', d_type);
    end
    if ~ismember(aggregation, VALID_AGGREGATIONS)
        error('build_cost_function:BadAggregation', ...
            'Unknown aggregation ''%s''. Valid: mean, max, min.', aggregation);
    end

    mask_data(k).lin_idx     = lin_idx;
    mask_data(k).sin_w       = sin_theta_ext(ext_t);
    mask_data(k).aggregation = aggregation;
    mask_data(k).d_type      = d_type;
    mask_data(k).d_weight    = d_weight;
end

cost_fn = @cost_fn_impl;

    function J = cost_fn_impl(x)
        % Decode x into complex weights according to the optimization mode.
        if strcmp(mode, 'amplitude_only')
            weights_complex = complex(x(:));         % real amplitudes, phase 0
        else
            weights_complex = x_to_weights(x);
            if strcmp(mode, 'phase_only')
                amplitudes = abs(weights_complex);
                safe = amplitudes;
                safe(~(amplitudes > 0)) = 1.0;       % avoid div-by-zero
                weights_complex = weights_complex ./ safe;
            end
        end

        % Power-normalise weights so cost is independent of global amplitude scale.
        w_power = sum(abs(weights_complex) .^ 2);
        weights_complex = weights_complex ./ sqrt(max(w_power, 1e-30));

        % Array factor (and total power for the secondary-stack case).
        af_primary = compute_array_factor(weights_complex, element_patterns_stacked);
        if ~isempty(element_patterns_secondary)
            af_secondary = compute_array_factor(weights_complex, element_patterns_secondary);
            power_grid = abs(af_primary) .^ 2 + abs(af_secondary) .^ 2;
        else
            power_grid = abs(af_primary) .^ 2;
        end

        % Composite cost: weighted sum over directives, sampled at sparse points.
        total_cost = 0.0;
        for kk = 1:numel(mask_data)
            masked_power = power_grid(mask_data(kk).lin_idx);   % (K x 1)
            term = directive_cost(masked_power, mask_data(kk).sin_w, ...
                                  mask_data(kk).d_type, mask_data(kk).aggregation);
            total_cost = total_cost + mask_data(kk).d_weight * term;
        end
        J = double(total_cost);
    end
end


function c = directive_cost(masked_power, sin_weights, directive_type, aggregation)
% Scalar cost contribution for one directive. Peak: -metric; Null: +metric.
% Ported from src/cost/cost_function.py:_directive_cost.
if isempty(masked_power)
    c = 0.0;
    return;
end

switch aggregation
    case 'mean'
        total_w = sum(sin_weights);
        metric  = dot(masked_power, sin_weights) / max(total_w, 1e-30);
    case 'max'
        metric  = max(masked_power);
    otherwise   % 'min'
        metric  = min(masked_power);
end

if strcmp(directive_type, 'peak')
    c = -metric;
else
    c = metric;
end
end
