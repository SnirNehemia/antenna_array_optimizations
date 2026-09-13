function summary = compare_approaches(results, scenario, array, settle_steps, results_dir)
% ══════════════════════════════════════════════════════════════════
% COMPARE_APPROACHES
% Score both approaches against what this array can actually achieve.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   summary = COMPARE_APPROACHES(results, scenario, array, settle_steps, ...
%                                results_dir)
%
%   THE SCORE.
%
%       score = 100 * (fraction of settled steps within 3 dB of the oracle)
%
%   It is measured against the ORACLE rather than against an absolute SINR
%   target, which means every array is scored against its own ceiling. That
%   matters because the arrays here differ by 20 dB in what they can achieve,
%   and an absolute target would simply rank them by aperture.
%
%   The alternative -- availability, the fraction of time above a fixed SINR --
%   saturates on easy geometries and cannot rank algorithms at all: the previous
%   work measured it moving 0.05 percentage points for a change worth 35 points
%   on this metric. It is not used here, not even as a secondary figure.
%
%   THE RULE THIS FILE EXISTS TO ENFORCE: the achievable dB is printed on the
%   same line as the score, always. A score of 99 on a two-element array whose
%   ceiling is 11 dB is not better than 72 on a sixteen-element array whose
%   ceiling is 31 dB, and a table of bare scores misleads everyone who reads
%   one. If you ever find yourself quoting a score without its ceiling, the
%   number has stopped meaning anything.
%
%   A pleasant consequence worth keeping: an array that can do nothing scores
%   100 for doing nothing, because its quiescent beam IS its optimum. The
%   one-element array is not marked as failing -- it is marked as achieving
%   everything available to it, which is the truth.
%
%   WHY SETTLED STEPS ONLY. The covariance needs a few horizons to fill and the
%   motion classifier needs its window. Scoring the warm-up measures the
%   estimator's convergence, not the algorithm, and mixing the two is how a run
%   length silently becomes a confound.
%
%   Inputs:
%       results      : struct from run_closed_loop.
%       scenario     : struct from make_scenario. TRUTH -- reporting only.
%       array        : struct from make_array.
%       settle_steps : first step to score from. Units: count.
%       results_dir  : folder for the CSV. Created if absent.
%
%   Outputs:
%       summary : struct with fields
%           names          : {1 x n} approach names.
%           mean_sinr_db   : (1 x n) mean settled SINR. Units: dB.
%           score          : (1 x n) percent of settled steps within 3 dB.
%           achievable_db  : scalar, the oracle's mean settled SINR. Units: dB.
%           gap_db         : (1 x n) achievable_db - mean_sinr_db. Units: dB.

% ────────────────────────── CONSTANTS ─────────────────────────────

% How close to the oracle counts as "keeping up". 3 dB is a factor of two in
% power -- the conventional engineering threshold for a difference that matters,
% and it is a stated judgement rather than a derivation.
TRACKING_TOLERANCE_DB = 3.0;

scored_steps = settle_steps : scenario.n_steps;
if isempty(scored_steps)
    error('compare_approaches:NothingToScore', ...
        ['settle_steps is %d but the scenario has only %d steps, so there is ' ...
         'nothing left to score.'], settle_steps, scenario.n_steps);
end

n_approaches  = numel(results.names);
oracle_scored = results.oracle_sinr_db(scored_steps);
achievable_db = mean(oracle_scored);

mean_sinr_db = zeros(1, n_approaches);
score        = zeros(1, n_approaches);
for i = 1:n_approaches
    approach_scored = results.sinr_db(scored_steps, i);
    mean_sinr_db(i) = mean(approach_scored);
    score(i)        = 100 * mean((oracle_scored - approach_scored) <= TRACKING_TOLERANCE_DB);
end

summary = struct( ...
    'names',         {results.names}, ...
    'mean_sinr_db',  mean_sinr_db, ...
    'score',         score, ...
    'achievable_db', achievable_db, ...
    'gap_db',        achievable_db - mean_sinr_db);

% ────────────────────────── PRINT ─────────────────────────────────

fprintf('\n%s  |  %s  |  %s\n', array.name, array.component, scenario.description);
fprintf('scored over steps %d-%d of %d; achievable (oracle) %.1f dB\n', ...
    scored_steps(1), scored_steps(end), scenario.n_steps, achievable_db);
fprintf('%-18s %12s %10s %12s %14s\n', ...
    'approach', 'mean SINR', 'score', 'gap to max', 'of achievable');
for i = 1:n_approaches
    fprintf('%-18s %9.1f dB %9.0f %9.1f dB %11.1f dB\n', ...
        results.names{i}, mean_sinr_db(i), score(i), summary.gap_db(i), achievable_db);
end

% Direction finding. The beamscan places the null; MUSIC is computed alongside
% it and reported, never used. On a moving jammer MUSIC sticks and jumps,
% because a source smeared across the covariance memory is not a point source --
% see estimate_jammer_angle for the mechanism.
true_theta_scored = scenario.jammer_theta_deg(scored_steps);
beamscan_error    = results.beamscan_theta_deg(scored_steps) - true_theta_scored;
music_error       = results.music_theta_deg(scored_steps)    - true_theta_scored;

fprintf('%-18s %12s %12s %12s\n', 'direction finding', 'median err', 'mean |err|', 'max |err|');
fprintf('%-18s %8.2f deg %8.2f deg %8.2f deg\n', '  beamscan (used)', ...
    median(beamscan_error, 'omitnan'), mean(abs(beamscan_error), 'omitnan'), ...
    max(abs(beamscan_error)));
fprintf('%-18s %8.2f deg %8.2f deg %8.2f deg\n', '  MUSIC (reported)', ...
    median(music_error, 'omitnan'), mean(abs(music_error), 'omitnan'), ...
    max(abs(music_error)));

% Approach 1's classification, which is reporting rather than scoring.
observed_behaviours = unique(results.behaviour(scored_steps));
fprintf('Approach 1 classified the jammer as: %s (truth: %s)\n', ...
    strjoin(observed_behaviours(:).', ', '), scenario.behaviour);

% ────────────────────────── SAVE ──────────────────────────────────

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

csv_path = fullfile(results_dir, sprintf('run_%s_%s.csv', ...
    matlab.lang.makeValidName(array.name), scenario.behaviour));

file_id = fopen(csv_path, 'w');
if file_id < 0
    error('compare_approaches:CannotWrite', 'Could not open %s for writing.', csv_path);
end

fprintf(file_id, ['step,true_jammer_theta_deg,jammer_on,beamscan_theta_deg,' ...
                  'music_theta_deg,aimed_theta_deg,behaviour']);
fprintf(file_id, ',%s', results.names{:});
fprintf(file_id, ',oracle\n');

for step_index = 1:scenario.n_steps
    fprintf(file_id, '%d,%.1f,%d,%.2f,%.2f,%.2f,%s', step_index, ...
        scenario.jammer_theta_deg(step_index), scenario.jammer_on(step_index), ...
        results.beamscan_theta_deg(step_index), results.music_theta_deg(step_index), ...
        results.jammer_theta_deg(step_index), results.behaviour{step_index});
    fprintf(file_id, ',%.2f', results.sinr_db(step_index, :));
    fprintf(file_id, ',%.2f\n', results.oracle_sinr_db(step_index));
end
fclose(file_id);

fprintf('saved %s\n', csv_path);
end
