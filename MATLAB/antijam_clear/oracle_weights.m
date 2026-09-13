function weights = oracle_weights(scenario, step_index, array)
% ══════════════════════════════════════════════════════════════════
% ORACLE_WEIGHTS
% The best any beamformer could do here. NOT DELIVERABLE -- it is told the truth.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   weights = ORACLE_WEIGHTS(scenario, step_index, array)
%
%   The same minimum-variance solve as max_sinr_weights, but built from an
%   EXACT covariance assembled from the true directions and the true powers:
%
%       R_true = sigma_s^2 e_s e_s' + sigma_j^2 e_j e_j' + sigma_n^2 I
%       w      = R_true \ e_s / (e_s' * (R_true \ e_s))
%
%   No estimation, no finite snapshots, no steering error, and therefore no
%   diagonal loading -- with the true steering vector there is no mismatch for
%   loading to protect against, and adding any would only make it worse.
%
%   WHY THIS FILE EXISTS. So that every score can be quoted as "82, against an
%   achievable 31.4 dB" rather than "82 out of 100". Without it, a 2-element
%   array scoring 99 because its ceiling is 11 dB looks better than a 16-element
%   array scoring 72 against a ceiling of 31 dB, which is exactly backwards. A
%   table of bare closeness scores misleads every reader who sees one.
%
%   It also gives the metric a pleasant property worth keeping: an array that
%   can do nothing scores 100 for doing nothing, because its quiescent beam IS
%   its optimum. The one-element array is not marked as failing; it is marked as
%   achieving everything available to it, which is the truth.
%
%   THIS IS NOT AN ALGORITHM. It takes the scenario struct -- the truth -- as
%   its first argument, which no algorithm in this folder may do. It is kept in
%   its own file so that is impossible to miss, and it exists only to bound the
%   others.
%
%   Inputs:
%       scenario   : struct from make_scenario. TRUTH.
%       step_index : 1-based step to evaluate. Units: count.
%       array      : struct from make_array.
%
%   Outputs:
%       weights : (n_elements x 1) complex, unit gain on the true signal
%                 direction. Units: dimensionless.

% ────────────────────────── TRUE STEERING VECTORS ─────────────────

signal_steering = steering_vector(array.element_patterns, array.theta_deg, ...
                                  array.phi_deg, scenario.signal_theta_deg, ...
                                  scenario.signal_phi_deg);

% ────────────────────────── THE EXACT COVARIANCE ──────────────────

% Each source contributes a rank-one term: its power times the outer product of
% its steering vector. This is what sample_covariance is estimating, written
% down exactly.
n_elements = array.n_elements;

covariance = scenario.signal_power * (signal_steering * signal_steering') ...
             + scenario.noise_power * eye(n_elements);

if scenario.jammer_on(step_index)
    jammer_steering = steering_vector(array.element_patterns, array.theta_deg, ...
                                      array.phi_deg, ...
                                      scenario.jammer_theta_deg(step_index), ...
                                      scenario.jammer_phi_deg(step_index));

    covariance = covariance ...
                 + scenario.jammer_power * (jammer_steering * jammer_steering');
end

% ────────────────────────── THE SOLVE ─────────────────────────────

whitened_steering = covariance \ signal_steering;
normalisation     = signal_steering' * whitened_steering;

weights = whitened_steering / conj(normalisation);

if ~all(isfinite(weights))
    error('oracle_weights:NonFiniteWeights', ...
        ['The oracle solve produced non-finite weights, which should be ' ...
         'impossible: the exact covariance is positive definite because it ' ...
         'contains sigma_n^2 * I. Check the scenario power levels.']);
end
end
