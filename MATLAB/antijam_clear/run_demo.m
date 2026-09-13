% ══════════════════════════════════════════════════════════════════
% RUN_DEMO
% Anti-jam beamforming, two approaches, read top to bottom.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
% THE PROBLEM. One jammer, direction unknown, somewhere on the sky. We want to
% keep receiving a wanted signal from a known direction. We get to choose one
% complex weight per antenna element, once per step.
%
% TWO WAYS TO DO IT, deliberately kept separate:
%
%   Approach 1, DETECT AND NULL. Work out where the jammer is and what it is
%   doing, then put a zero in the radiation pattern there. Explicit and
%   geometric: it needs two directions and no data at all.
%
%   Approach 2, MAX SINR. Never ask where the jammer is. Solve directly for the
%   weights that let through as little non-signal power as possible. Implicit,
%   closed form, one line of algebra.
%
% Run it: run_demo
%
% ══════════════════════════════════════════════════════════════════

clear; clc;

this_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(this_folder, '..', 'matlab_utils'));
addpath(this_folder);

results_dir = fullfile(this_folder, 'results');


%% ────────────────────── THE CONSTANTS ────────────────────────────
%
% Two tuned numbers -- FORGETTING_LAMBDA and LOADING_FACTOR -- and everything
% else derived from them. They all fit on this screen.

% How far back the covariance remembers. Both memory numbers below follow from
% it, so if this moves they move with it.
FORGETTING_LAMBDA = 0.90;

% The covariance weights a block from k steps ago by (1-lambda)*lambda^k. That
% decaying weighting has TWO useful one-number summaries, and they are not the
% same number -- the distinction is worth about 1 dB on a drifting jammer.
%
% How much smoothing this is equivalent to: 10 steps. Sets how long the array
% remembers, why an on/off edge cannot be resolved from this covariance, how
% long a null survives the jammer switching off, and the motion window.
COVARIANCE_HORIZON_STEPS = round(1 / (1 - FORGETTING_LAMBDA));

% The MEAN AGE of the data in that weighted average: 9 steps. The beamscan
% reports where the jammer was on average, and that average is this old -- so
% this, not the window length above, is the lag, and this is what the null is
% led by. Measured directly: at 0.10 and 0.20 deg/step the lag divided by the
% rate comes out at 9.00 steps exactly.
ESTIMATE_LAG_STEPS = round(FORGETTING_LAMBDA / (1 - FORGETTING_LAMBDA));

% Diagonal loading, as a multiple of the estimated noise floor. Chosen by
% measuring SINR against STEERING MISMATCH rather than against null depth: see
% results/loading_sweep.csv. 30 is the best worst case over 0-3 degrees of
% mismatch and 0-20 dB of signal strength, and costs at most 1.6 dB against
% tuning it separately for each condition.
LOADING_FACTOR = 30;

% Skip the warm-up when scoring: the covariance needs a few horizons to fill and
% the motion classifier needs its window. Scoring the warm-up measures how fast
% the estimator converges, not how well the algorithm works.
SETTLE_STEPS = 3 * COVARIANCE_HORIZON_STEPS;


%% ────────────────────── LOAD THE ARRAY ───────────────────────────
%
% Real CST element patterns. The component must be named: the patch arrays
% export Copol/Cross, the dipole and monopole arrays export Theta/Phi, and there
% is no default that is right for both.

array_folder     = fullfile(this_folder, '..', '..', 'data', 'spacing0.6');
component_name   = 'Copol';
signal_theta_deg = 30;
signal_phi_deg   = 0;

array = make_array(array_folder, component_name, signal_theta_deg, signal_phi_deg);

% Everything below depends on these measured facts, so print them first.
fprintf('ARRAY: %s (%s)\n', array.name, array.component);
fprintf('  %d elements, %d degrees of freedom after the signal constraint\n', ...
    array.profile.n_elements, array.profile.degrees_of_freedom);
fprintf('  quiescent directivity at the target : %.1f dBi\n', array.profile.quiescent_dbi);
fprintf('  3 dB beamwidth                      : %.1f deg (theta), %.1f deg (phi)\n', ...
    array.profile.hpbw_theta_deg, array.profile.hpbw_phi_deg);
fprintf('  GUARD SECTOR (derived, half the wider beamwidth) : %.1f deg\n', ...
    array.profile.guard_deg);
fprintf('  mirror ambiguity e(theta) vs e(180-theta)        : %.4f%s\n', ...
    array.profile.mirror_coherence, ambiguity_note(array.profile.is_mirror_ambiguous));
fprintf('  snapshots per step (2 per element)               : %d\n', ...
    array.snapshots_per_step);
fprintf('  covariance memory: %d steps equivalent window, %d steps mean age\n\n', ...
    COVARIANCE_HORIZON_STEPS, ESTIMATE_LAG_STEPS);


%% ────────────────────── ONE INSTANT, STEP BY STEP ────────────────
%
% The whole method at a single moment, with nothing hidden in a loop. This is
% the part to read aloud.

scenario = make_scenario(array.profile, 'steady', signal_theta_deg, signal_phi_deg, ...
                         25, ...      % jammer separation [deg] - must exceed the guard
                         0, ...       % signal-to-noise [dB]
                         20);         % jammer-to-noise [dB]

fprintf('SCENARIO: %s\n', scenario.description);
fprintf('  the true jammer angle lives in this struct and is never passed to an algorithm\n\n');

% The receiver collects K complex samples from every element.
snapshots = simulate_snapshots(scenario, 1, array);

% What the array has been hearing. On the first step there is no past to blend.
covariance = sample_covariance([], snapshots, FORGETTING_LAMBDA);

% APPROACH 1, part one: where is the jammer, and what is it doing?
jammer_state = detect_jammer(covariance, [], array, COVARIANCE_HORIZON_STEPS, ...
                             ESTIMATE_LAG_STEPS);

fprintf('DETECTED: jammer %s at theta = %.1f deg, classified %s\n', ...
    presence_word(jammer_state.is_present), jammer_state.theta_deg, jammer_state.behaviour);
fprintf('  (truth, for our eyes only: %.1f deg)\n', scenario.jammer_theta_deg(1));
fprintf('  %s\n\n', jammer_state.reason);

% APPROACH 1, part two: hold the signal, put a zero on the jammer.
[weights_detect_and_null, null_reason] = detect_and_null(array, jammer_state);

% APPROACH 2: never look for the jammer; just refuse everything that is not signal.
weights_max_sinr = max_sinr_weights(covariance, array, LOADING_FACTOR);

fprintf('APPROACH 1 (detect and null): %s\n', null_reason);
fprintf('  SINR %.1f dB\n', output_sinr_db(weights_detect_and_null, scenario, 1, array));
fprintf('APPROACH 2 (max SINR): no jammer angle used; diagonal loading %g x noise floor\n', ...
    LOADING_FACTOR);
fprintf('  SINR %.1f dB\n', output_sinr_db(weights_max_sinr, scenario, 1, array));
fprintf('ORACLE (told the truth, not deliverable): SINR %.1f dB\n\n', ...
    output_sinr_db(oracle_weights(scenario, 1, array), scenario, 1, array));


%% ────────────────────── THE WHOLE RUN, THREE BEHAVIOURS ──────────
%
% Exactly the five calls above, once per step. run_closed_loop carries only two
% things between steps: the covariance, and the detector's history of angles.

behaviours = {'steady', 'onoff', 'drift'};

for behaviour_index = 1:numel(behaviours)

    scenario = make_scenario(array.profile, behaviours{behaviour_index}, ...
                             signal_theta_deg, signal_phi_deg, 25, 0, 20);

    results = run_closed_loop(array, scenario, FORGETTING_LAMBDA, ...
                              COVARIANCE_HORIZON_STEPS, ESTIMATE_LAG_STEPS, ...
                              LOADING_FACTOR);

    compare_approaches(results, scenario, array, SETTLE_STEPS, results_dir);
    plot_pattern_cut(results, array, scenario, scenario.n_steps, results_dir);
end


%% ────────────────────── ARRAYS THAT CANNOT DO THIS ───────────────
%
% Two of the arrays in data/ provably cannot null. The method should say so
% clearly rather than crash or emit a meaningless number - knowing its own
% limits is part of what makes it trustworthy.

fprintf('\n\nARRAYS AT THEIR LIMITS\n');

limited_folders    = {'Dipole', 'patch_back2back'};
limited_components = {'Theta',  'Copol'};

for limited_index = 1:numel(limited_folders)

    limited_array = make_array(fullfile(this_folder, '..', '..', 'data', ...
                                        limited_folders{limited_index}), ...
                               limited_components{limited_index}, ...
                               signal_theta_deg, signal_phi_deg);

    fprintf('\n%s: %d element(s), %d degrees of freedom\n', ...
        limited_folders{limited_index}, limited_array.profile.n_elements, ...
        limited_array.profile.degrees_of_freedom);

    if ~limited_array.profile.is_target_illuminated
        fprintf('  the array is %.1f dB down at the target: that is a pattern null, not a\n', ...
            -limited_array.profile.target_visibility_db);
        fprintf('  beam. It cannot receive its own signal there, so this is not a test.\n');
        continue
    end

    limited_state = detect_jammer(eye(limited_array.n_elements), [], limited_array, ...
                                  COVARIANCE_HORIZON_STEPS, ESTIMATE_LAG_STEPS);
    [~, limited_reason] = detect_and_null(limited_array, limited_state);
    fprintf('  %s\n', limited_reason);
end

fprintf('\nDone. Results and figures are in %s\n', results_dir);


%% ────────────────────── SMALL PRINTING HELPERS ───────────────────

function note = ambiguity_note(is_ambiguous)
% Spell out what the mirror coherence number means, so the printout is readable.
if is_ambiguous
    note = '  <- AMBIGUOUS: harmless for nulling, must be folded before tracking';
else
    note = '  (unambiguous)';
end
end


function word = presence_word(is_present)
% 'present' / 'not present', for the printout.
if is_present
    word = 'present';
else
    word = 'NOT present';
end
end
