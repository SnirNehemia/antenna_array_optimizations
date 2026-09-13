function results = run_closed_loop(array, scenario, forgetting_lambda, ...
                                   covariance_horizon_steps, estimate_lag_steps, ...
                                   loading_factor)
% ══════════════════════════════════════════════════════════════════
% RUN_CLOSED_LOOP
% The whole method, once per time step. Five named calls and nothing else.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   results = RUN_CLOSED_LOOP(array, scenario, forgetting_lambda, ...
%                             covariance_horizon_steps, estimate_lag_steps, ...
%                             loading_factor)
%
%   Steps through a scenario, and at each step:
%
%       1. simulate_snapshots  -- the receiver collects K samples per element
%       2. sample_covariance   -- update what the array has been hearing lately
%       3. detect_jammer       -- Approach 1, part one: where is it, what is it doing
%       4. detect_and_null     -- Approach 1, part two: hold signal, zero the jammer
%       5. max_sinr_weights    -- Approach 2: minimise everything that is not signal
%
%   plus the quiescent beam as a baseline and the oracle as the ceiling, both
%   scored by output_sinr_db. That is the entire method; everything else in this
%   folder is either feeding these five calls or reporting on them.
%
%   WHAT IS CARRIED BETWEEN STEPS, and nothing else is: the covariance, and the
%   detector's history of past angle estimates. Both are returned by the
%   function that owns them and passed back in, so the state is visible in the
%   loop rather than hidden in a persistent variable.
%
%   THE TRUTH BOUNDARY. scenario is passed to exactly three things here:
%   simulate_snapshots (to build the data), oracle_weights (which is the bound,
%   not a deliverable) and output_sinr_db (which scores). It is never passed to
%   detect_jammer, detect_and_null or max_sinr_weights -- check the three calls
%   below and you will see they receive only the covariance and the array.
%
%   Inputs:
%       array                    : struct from make_array.
%       scenario                 : struct from make_scenario. TRUTH.
%       forgetting_lambda        : covariance memory factor. Units: dimensionless.
%       covariance_horizon_steps : 1/(1 - forgetting_lambda), the effective
%                                  window length. Units: steps.
%       estimate_lag_steps       : forgetting_lambda/(1 - forgetting_lambda),
%                                  the mean age of the covariance. Units: steps.
%       loading_factor           : diagonal loading, as a multiple of the noise
%                                  floor. Units: dimensionless.
%
%   Outputs:
%       results : struct with fields
%           names          : {1 x 3} cell, the approach names.
%           sinr_db        : (n_steps x 3) SINR per step per approach. Units: dB.
%           oracle_sinr_db : (n_steps x 1) the achievable ceiling. Units: dB.
%           jammer_theta_deg : (n_steps x 1) where Approach 1 aimed its null.
%                              Units: degrees.
%           beamscan_theta_deg : (n_steps x 1) raw beamscan estimate, which is
%                              what places the null. Units: degrees.
%           music_theta_deg  : (n_steps x 1) raw MUSIC estimate, reported for
%                              comparison only. Units: degrees.
%           behaviour      : {n_steps x 1} cell, Approach 1's classification.
%           final_weights  : {1 x 3} cell, the weights at the last step, for
%                            plotting.

n_steps = scenario.n_steps;

results = struct( ...
    'names',            {{'quiescent', 'detect_and_null', 'max_sinr'}}, ...
    'sinr_db',          nan(n_steps, 3), ...
    'oracle_sinr_db',   nan(n_steps, 1), ...
    'jammer_theta_deg', nan(n_steps, 1), ...
    'beamscan_theta_deg', nan(n_steps, 1), ...
    'music_theta_deg',  nan(n_steps, 1), ...
    'behaviour',        {cell(n_steps, 1)}, ...
    'final_weights',    {cell(1, 3)});

% The quiescent beam never changes: it does not adapt, which is the point of it.
signal_steering  = steering_vector(array.element_patterns, array.theta_deg, ...
                                   array.phi_deg, array.profile.signal_theta_deg, ...
                                   array.profile.signal_phi_deg);
quiescent_beam   = quiescent_weights(signal_steering);

% ────────────────────────── THE LOOP ──────────────────────────────

covariance       = [];      % grown by sample_covariance
detector_history = [];      % grown by detect_jammer

for step_index = 1:n_steps

    % 1. What the receiver hears this step.
    snapshots  = simulate_snapshots(scenario, step_index, array);

    % 2. What it has been hearing lately.
    covariance = sample_covariance(covariance, snapshots, forgetting_lambda);

    % 3-4. Approach 1: find the jammer, then null it.
    [jammer_state, detector_history] = detect_jammer(covariance, detector_history, ...
                                                     array, covariance_horizon_steps, ...
                                                     estimate_lag_steps);
    weights_detect_and_null = detect_and_null(array, jammer_state);

    % 5. Approach 2: never look for the jammer at all.
    weights_max_sinr = max_sinr_weights(covariance, array, loading_factor);

    % Scoring. These see the truth; the three calls above do not.
    results.sinr_db(step_index, 1) = output_sinr_db(quiescent_beam, scenario, step_index, array);
    results.sinr_db(step_index, 2) = output_sinr_db(weights_detect_and_null, scenario, step_index, array);
    results.sinr_db(step_index, 3) = output_sinr_db(weights_max_sinr, scenario, step_index, array);
    results.oracle_sinr_db(step_index) = output_sinr_db(oracle_weights(scenario, step_index, array), ...
                                                        scenario, step_index, array);

    results.jammer_theta_deg(step_index)   = jammer_state.theta_deg;
    results.beamscan_theta_deg(step_index) = jammer_state.measured_theta_deg;
    results.music_theta_deg(step_index)    = jammer_state.music_theta_deg;
    results.behaviour{step_index}        = jammer_state.behaviour;

    if step_index == n_steps
        results.final_weights = {quiescent_beam, weights_detect_and_null, weights_max_sinr};
    end
end

% A run that completes with non-finite scores looks like a good run and is not
% one. Assert the values, not merely that the loop finished.
if ~all(isfinite(results.sinr_db(:))) || ~all(isfinite(results.oracle_sinr_db))
    error('run_closed_loop:NonFiniteScores', ...
        ['The run completed but produced non-finite SINR values. This is the ' ...
         'failure mode that looks like success -- do not ignore it.']);
end
end
