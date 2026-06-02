function result = run_optimizer(element_patterns_stacked, theta_deg, phi_deg, ...
                                directives, optimizer_config, element_patterns_secondary)
% RUN_OPTIMIZER  Multi-start bound-constrained optimization of complex weights.
%
%   result = RUN_OPTIMIZER(element_patterns_stacked, theta_deg, phi_deg, ...
%               directives, optimizer_config, element_patterns_secondary)
%
%   MATLAB port of the scipy L-BFGS-B optimizer using fmincon (Optimization
%   Toolbox). Runs n_restarts user restarts (the first uniform when
%   use_uniform_init) followed by one single-element restart per element when
%   use_single_element_init. Returns the lowest-cost solution found.
%
%   NOTE (determinism): fmincon is not scipy's L-BFGS-B, and MATLAB rng(0)
%   (Mersenne Twister) is not numpy default_rng(0) (PCG64), so iterates and
%   random restarts do NOT match Python step-for-step. Both minimise the same
%   cost J(x), so the converged cost/pattern agree on well-posed problems.
%
%   Inputs:
%       element_patterns_stacked   : (N_elements, N_theta, N_phi) complex [V/m].
%       theta_deg, phi_deg         : angle grids [deg].
%       directives                 : cell array of directive structs.
%       optimizer_config           : struct. Required: max_iterations,
%                                    cost_tolerance, n_restarts. Optional:
%                                    amplitude_bounds ([] or [a_min a_max]),
%                                    phase_only, amplitude_only, use_uniform_init,
%                                    use_single_element_init.
%       element_patterns_secondary : optional cross-pol stack for total-power. [].
%
%   Outputs:
%       result : struct with fields weights_complex (N x 1), cost_history (vector),
%                result (struct: x, fun, success), all_cost_histories (cell),
%                all_run_labels (cell), best_run_index (1-BASED), weights_history (cell).
%
%   Ported from src/optimize/optimizer.py:run_optimizer.
%   Part of: Antenna Array Pattern Optimization Tool.

if nargin < 6
    element_patterns_secondary = [];
end

REQUIRED = {'max_iterations', 'cost_tolerance', 'n_restarts'};
for i = 1:numel(REQUIRED)
    if ~isfield(optimizer_config, REQUIRED{i})
        error('run_optimizer:MissingKey', ...
            'Missing required optimizer config key: ''%s''.', REQUIRED{i});
    end
end

max_iterations = optimizer_config.max_iterations;
cost_tolerance = optimizer_config.cost_tolerance;
n_restarts     = optimizer_config.n_restarts;
amplitude_bounds        = cfg_get(optimizer_config, 'amplitude_bounds', []);
phase_only              = cfg_get(optimizer_config, 'phase_only', false);
amplitude_only          = cfg_get(optimizer_config, 'amplitude_only', false);
use_uniform_init        = cfg_get(optimizer_config, 'use_uniform_init', true);
use_single_element_init = cfg_get(optimizer_config, 'use_single_element_init', true);

if phase_only && amplitude_only
    error('run_optimizer:ModeConflict', ...
        '''phase_only'' and ''amplitude_only'' cannot both be true.');
end
if phase_only
    mode = 'phase_only';
elseif amplitude_only
    mode = 'amplitude_only';
else
    mode = 'standard';
end

n_elements = size(element_patterns_stacked, 1);

cost_fn   = build_cost_function(element_patterns_stacked, theta_deg, phi_deg, ...
                                directives, mode, element_patterns_secondary);
[lb, ub]  = build_bounds(n_elements, amplitude_bounds, mode);

% Seeded RNG so restarts are reproducible across MATLAB runs.
rng(0);

best_result        = [];
best_cost_history  = [];
best_xk_history    = {};
best_final_cost    = inf;
best_run_index     = 0;       % 1-based once set
all_cost_histories = {};
all_run_labels     = {};
run_index          = 0;

% ── Phase 1: user-configured restarts ─────────────────────────────
for restart_index = 0:(n_restarts - 1)
    if restart_index == 0 && use_uniform_init
        x_initial = uniform_initial_x(n_elements, mode);
        run_label = 'Uniform init';
    else
        x_initial = random_initial_x(n_elements, mode);
        run_label = sprintf('Random init %d', restart_index);
    end
    run_index = run_index + 1;

    [xopt, fval, success, cost_history, xk_history] = run_single( ...
        cost_fn, x_initial, lb, ub, max_iterations, cost_tolerance);
    all_cost_histories{end + 1, 1} = cost_history; %#ok<AGROW>
    all_run_labels{end + 1, 1}     = run_label;    %#ok<AGROW>

    if fval < best_final_cost
        best_final_cost   = fval;
        best_result       = struct('x', xopt, 'fun', fval, 'success', success);
        best_cost_history = cost_history;
        best_xk_history   = xk_history;
        best_run_index    = run_index;
    end
end

% ── Phase 2: single-element initializations ───────────────────────
if use_single_element_init
    for element_idx = 0:(n_elements - 1)
        x_initial = single_element_initial_x(n_elements, element_idx, mode);
        run_label = sprintf('Element %d init', element_idx);
        run_index = run_index + 1;

        [xopt, fval, success, cost_history, xk_history] = run_single( ...
            cost_fn, x_initial, lb, ub, max_iterations, cost_tolerance);
        all_cost_histories{end + 1, 1} = cost_history; %#ok<AGROW>
        all_run_labels{end + 1, 1}     = run_label;    %#ok<AGROW>

        if fval < best_final_cost
            best_final_cost   = fval;
            best_result       = struct('x', xopt, 'fun', fval, 'success', success);
            best_cost_history = cost_history;
            best_xk_history   = xk_history;
            best_run_index    = run_index;
        end
    end
end

if isempty(best_result)
    error('run_optimizer:NoRuns', ...
        'No optimization runs were performed. Set n_restarts >= 1 or enable single-element init.');
end

% ── Decode best solution (no power-normalisation here) ────────────
weights_complex = decode_x(best_result.x, mode);

weights_history = cell(numel(best_xk_history), 1);
for i = 1:numel(best_xk_history)
    weights_history{i} = decode_x(best_xk_history{i}, mode);
end

result = struct();
result.weights_complex    = weights_complex;
result.cost_history       = best_cost_history;
result.result             = best_result;
result.all_cost_histories = all_cost_histories;
result.all_run_labels     = all_run_labels;
result.best_run_index     = best_run_index;     % 1-based (Python is 0-based)
result.weights_history    = weights_history;
end


% ────────────────────────── HELPERS ───────────────────────────────

function val = cfg_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end


function [lb, ub] = build_bounds(n_elements, amplitude_bounds, mode)
% Ported from src/optimize/optimizer.py:_build_lbfgsb_bounds.
if strcmp(mode, 'phase_only') || isempty(amplitude_bounds)
    lb = [];
    ub = [];
    return;
end
a_min = amplitude_bounds(1);
a_max = amplitude_bounds(2);
if strcmp(mode, 'amplitude_only')
    lb = a_min * ones(n_elements, 1);
    ub = a_max * ones(n_elements, 1);
else
    % Standard: bound each Re/Im component to [-a_max, a_max].
    lb = -a_max * ones(2 * n_elements, 1);
    ub =  a_max * ones(2 * n_elements, 1);
end
end


function x0 = uniform_initial_x(n_elements, mode)
if strcmp(mode, 'amplitude_only')
    x0 = ones(n_elements, 1);
else
    x0 = weights_to_x(ones(n_elements, 1));
end
end


function x0 = random_initial_x(n_elements, mode)
if strcmp(mode, 'amplitude_only')
    x0 = rand(n_elements, 1);
else
    random_phases_rad = rand(n_elements, 1) * 2 * pi;
    x0 = weights_to_x(exp(1i * random_phases_rad));
end
end


function x0 = single_element_initial_x(n_elements, element_idx0, mode)
% element_idx0 is 0-based (mirrors Python); +1 for MATLAB indexing.
if strcmp(mode, 'amplitude_only')
    x0 = zeros(n_elements, 1);
    x0(element_idx0 + 1) = 1.0;
else
    w = zeros(n_elements, 1);
    w(element_idx0 + 1) = 1.0;
    x0 = weights_to_x(w);
end
end


function w = decode_x(x, mode)
% Decode raw variable vector to complex weights (no power-normalisation).
if strcmp(mode, 'amplitude_only')
    w = complex(x(:));
    return;
end
w = x_to_weights(x);
if strcmp(mode, 'phase_only')
    a = abs(w);
    safe = a;
    safe(~(a > 0)) = 1.0;
    w = w ./ safe;
end
end


function [xopt, fval, success, cost_history, xk_history] = run_single( ...
        cost_fn, x_initial, lb, ub, max_iterations, cost_tolerance)
% Single fmincon run with per-iteration cost/x recording via OutputFcn.
cost_history = [];
xk_history   = {};

options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'MaxIterations', max_iterations, ...
    'MaxFunctionEvaluations', 1e6, ...
    'FunctionTolerance', cost_tolerance, ...
    'OptimalityTolerance', cost_tolerance, ...
    'StepTolerance', 1e-12, ...
    'OutputFcn', @record_iterate);

[xopt, fval, exitflag] = fmincon(cost_fn, x_initial, [], [], [], [], lb, ub, [], options);
success = exitflag > 0;
xopt    = xopt(:);

    function stop = record_iterate(x, optimValues, state)
        stop = false;
        if strcmp(state, 'iter')
            cost_history(end + 1, 1) = optimValues.fval;
            xk_history{end + 1, 1}   = x(:);
        end
    end
end
