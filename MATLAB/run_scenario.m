function results = run_scenario(scenario, element_patterns_complex, theta_deg, phi_deg, ...
                                n_side, d_over_lambda, cfg, element_patterns_secondary)
% RUN_SCENARIO  Classical and optimised patterns for one comparison scenario.
%
%   results = RUN_SCENARIO(scenario, element_patterns_complex, theta_deg, phi_deg, ...
%               n_side, d_over_lambda, cfg, element_patterns_secondary)
%
%   Computes, for each classical technique and the L-BFGS-B-style optimiser, the
%   power-normalised principal-plane cut and per-directive min/mean/max metrics.
%
%   Inputs:
%       scenario                   : struct with steer_theta_deg, steer_phi_deg,
%                                    directives (cell of directive structs).
%       element_patterns_complex   : (N_elements,N_theta,N_phi) complex [V/m].
%       theta_deg, phi_deg         : angle grids. Units: degrees.
%       n_side, d_over_lambda      : array geometry.
%       cfg                        : config struct (supplies optimizer settings).
%       element_patterns_secondary : optional cross-pol stack ([] -> copol only).
%
%   Outputs:
%       results : struct keyed by technique name (uniform, hamming, ..., optimized),
%                 each with fields weights_complex, cut_vm, directive_metrics.
%
%   Ported from scripts/compare_classical.py:run_scenario.
%   Part of: Antenna Array Pattern Optimization Tool.

if nargin < 8
    element_patterns_secondary = [];
end

CLASSICAL_TECHNIQUE_NAMES = {'uniform', 'hamming', 'hanning', 'kaiser_b3', ...
    'kaiser_b6', 'chebyshev_25dB', 'chebyshev_40dB', 'taylor_25dB'};

steer_theta_deg = scenario.steer_theta_deg;
steer_phi_deg   = scenario.steer_phi_deg;
directives      = scenario.directives;

results = struct();

% ── Classical techniques ──────────────────────────────────────────
for ti = 1:numel(CLASSICAL_TECHNIQUE_NAMES)
    name = CLASSICAL_TECHNIQUE_NAMES{ti};
    weights_raw = classical_weights(name, n_side, d_over_lambda, ...
        steer_theta_deg, steer_phi_deg, element_patterns_complex, theta_deg, phi_deg);
    [w, cut_vm, dm] = process_weights(weights_raw);
    results.(name) = struct('weights_complex', w, 'cut_vm', cut_vm, 'directive_metrics', dm);
end

% ── Optimised technique ───────────────────────────────────────────
n_restarts = cfg_get(cfg, 'n_restarts', 5);
optimizer_config = struct( ...
    'max_iterations',          cfg_get(cfg, 'max_iterations', 500), ...
    'cost_tolerance',          cfg_get(cfg, 'cost_tolerance', 1e-9), ...
    'n_restarts',              n_restarts, ...
    'use_uniform_init',        true, ...
    'use_single_element_init', false);
opt_result = run_optimizer(element_patterns_complex, theta_deg, phi_deg, ...
    directives, optimizer_config, element_patterns_secondary);
[w_opt, cut_opt, dm_opt] = process_weights(opt_result.weights_complex);
results.optimized = struct('weights_complex', w_opt, 'cut_vm', cut_opt, 'directive_metrics', dm_opt);


    function [weights_complex, cut_vm, dir_metrics] = process_weights(weights_complex)
        % Power-normalise, compute |AF| (or total-power magnitude), cut, metrics.
        weights_norm = power_normalize_weights(weights_complex);
        af_copol     = compute_array_factor(weights_norm, element_patterns_complex);
        if ~isempty(element_patterns_secondary)
            af_cross  = compute_array_factor(weights_norm, element_patterns_secondary);
            af_mag_vm = sqrt(abs(af_copol) .^ 2 + abs(af_cross) .^ 2);
        else
            af_mag_vm = abs(af_copol);
        end
        cut_vm = principal_plane_cut(af_mag_vm, phi_deg);
        dir_metrics = evaluate_directive_metrics(af_mag_vm, theta_deg, phi_deg, directives{1});
        for di = 2:numel(directives)
            dir_metrics(di) = evaluate_directive_metrics(af_mag_vm, theta_deg, phi_deg, directives{di});
        end
    end
end


function val = cfg_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end
