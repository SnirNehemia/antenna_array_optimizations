function sinr_db = output_sinr_db(weights, scenario, step_index, array)
% ══════════════════════════════════════════════════════════════════
% OUTPUT_SINR_DB
% The one scorer: what a given set of weights actually achieves.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   sinr_db = OUTPUT_SINR_DB(weights, scenario, step_index, array)
%
%   Everything in this folder is judged by this function and by nothing else,
%   so there is exactly one definition of "good" and no way for two parts of
%   the code to disagree about it.
%
%       P_signal = sigma_s^2 * |w' e_s|^2
%       P_jammer = sigma_j^2 * |w' e_j|^2       (zero while the jammer is off)
%       P_noise  = sigma_n^2 * ||w||^2
%       SINR     = P_signal / (P_jammer + P_noise)
%
%   Read it as three statements about the pattern the weights produce: how much
%   gain it has towards the signal, how much it has towards the jammer, and how
%   much noise it lets through -- the last being proportional to ||w||^2 because
%   receiver noise is independent on every element, so it adds in power rather
%   than in field.
%
%   This is a PATTERN-LEVEL SINR: it evaluates the weights against the true
%   directions and powers, rather than averaging over a finite draw of
%   snapshots. It is therefore the SINR those weights actually deliver, with no
%   estimation noise of its own -- which is what you want in a scorer, so that
%   run-to-run differences come from the algorithms and not from the scoring.
%
%   IT SEES THE TRUTH, so like oracle_weights it takes the scenario struct and
%   may never be called from inside an algorithm -- only from the comparison
%   and reporting code.
%
%   Inputs:
%       weights    : (n_elements x 1) complex weights to score.
%       scenario   : struct from make_scenario. TRUTH.
%       step_index : 1-based step to evaluate. Units: count.
%       array      : struct from make_array.
%
%   Outputs:
%       sinr_db : output signal-to-interference-plus-noise ratio. Units: dB.

weights = weights(:);

signal_steering = steering_vector(array.element_patterns, array.theta_deg, ...
                                  array.phi_deg, scenario.signal_theta_deg, ...
                                  scenario.signal_phi_deg);

% ────────────────────────── THE THREE POWERS ──────────────────────

signal_power = scenario.signal_power * abs(weights' * signal_steering) ^ 2;

jammer_power = 0.0;
if scenario.jammer_on(step_index)
    jammer_steering = steering_vector(array.element_patterns, array.theta_deg, ...
                                      array.phi_deg, ...
                                      scenario.jammer_theta_deg(step_index), ...
                                      scenario.jammer_phi_deg(step_index));

    jammer_power = scenario.jammer_power * abs(weights' * jammer_steering) ^ 2;
end

% Noise is independent per element, so it accumulates as the squared norm of
% the weights rather than coherently.
noise_power = scenario.noise_power * real(weights' * weights);

% ────────────────────────── THE RATIO ─────────────────────────────

interference_plus_noise = jammer_power + noise_power;
if interference_plus_noise <= 0
    error('output_sinr_db:NoDenominator', ...
        ['Interference plus noise came out as %g. Every run has receiver ' ...
         'noise, so this cannot legitimately be zero -- the weights are ' ...
         'probably all zero.'], interference_plus_noise);
end

sinr_db = 10 * log10(signal_power / interference_plus_noise);

if ~isfinite(sinr_db)
    error('output_sinr_db:NonFiniteSinr', ...
        ['SINR came out non-finite (signal power %g, interference plus noise ' ...
         '%g). A run that completes with a non-finite score looks like a good ' ...
         'run and is not one.'], signal_power, interference_plus_noise);
end
end
