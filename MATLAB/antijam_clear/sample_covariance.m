function covariance = sample_covariance(previous_covariance, snapshot_block, forgetting_lambda)
% ══════════════════════════════════════════════════════════════════
% SAMPLE_COVARIANCE
% What the array has been hearing lately, and how far back "lately" reaches.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   covariance = SAMPLE_COVARIANCE(previous_covariance, snapshot_block, ...
%                                  forgetting_lambda)
%
%   Two lines of real content:
%
%       R_block = (X * X') / K                       this step's own estimate
%       R       = lambda * R_previous + (1 - lambda) * R_block
%
%   R is an (n_elements x n_elements) matrix whose (m, n) entry is the average
%   correlation between what element m heard and what element n heard. A strong
%   source from one direction imprints a rank-one pattern on it, and that
%   imprint is everything both algorithms work from -- max_sinr_weights inverts
%   it, detect_jammer eigendecomposes it. Neither ever sees a raw snapshot again.
%
%   WHY FORGET AT ALL. The jammer moves and switches off. A plain running
%   average over all time would remember every place the jammer has ever been
%   and put the null at the average of them, which is nowhere useful. lambda
%   sets how far back the beamformer remembers.
%
%   WHY THE HORIZON IS THE NUMBER TO QUOTE. The weight given to a step k updates
%   ago falls as lambda^k, so the effective memory is
%
%       horizon = 1 / (1 - lambda)  steps
%
%   At lambda = 0.90 that is 10 steps, and this single number explains a
%   surprising share of the system's behaviour:
%
%       - a drifting jammer is nulled roughly where it was 10 steps ago, so the
%         null trails it by (drift rate) x 10 degrees;
%       - a detector reading THIS covariance cannot resolve an on/off edge,
%         because the covariance is smoothing the very transition being looked
%         for -- ten steps of smoothing across an edge that lasts one step;
%       - a null usefully survives the jammer switching off, for about that long.
%
%   Quote the horizon, not lambda. "The system remembers ten steps" is a
%   sentence anyone in the room can reason about.
%
%   THE IMPLEMENTATION DETAIL THAT IS LOAD-BEARING. The update is applied ONCE
%   per step, to the batch-averaged block. Applying it once per snapshot column
%   instead runs the recursion K times per step, shortening the true horizon by
%   a factor of K while lambda still reads 0.90 in the code. The previous
%   implementation had exactly that bug, and it made a lambda of 0.98 appear
%   correct. The form above is the correct one.
%
%   [FUTURE / DELIBERATELY OMITTED] The previous stack kept a SECOND covariance
%   at a much shorter memory (lambda = 0.5), used only for presence detection
%   and never to form weights, on the principle that a detector slower than the
%   beamformer cannot resolve edges the beamformer already smooths. It measurably
%   worked: presence error against the true duty cycle fell from 0.46 / 0.32 /
%   0.14 to 0.13 / 0.06 / 0.03 at toggle periods of 4 / 10 / 25 s. It is left out
%   here to keep one memory constant in the story, at the cost of on/off
%   detection being reactive and roughly one horizon late. Adding it back is
%   cheap and self-contained: call this function a second time with its own
%   previous_covariance and its own lambda, and feed only the presence test.
%   See README.md, "Deliberate omissions".
%
%   Inputs:
%       previous_covariance : (n_elements x n_elements) complex from the last
%                             step, or [] on the first step. Units: power.
%       snapshot_block      : (n_elements x K) complex from simulate_snapshots.
%       forgetting_lambda   : memory factor in [0, 1). Units: dimensionless.
%
%   Outputs:
%       covariance : (n_elements x n_elements) Hermitian positive definite.
%                    Units: power.

if ~(forgetting_lambda >= 0 && forgetting_lambda < 1)
    error('sample_covariance:BadLambda', ...
        ['forgetting_lambda must be in [0, 1); got %g. The memory horizon is ' ...
         '1/(1 - lambda) steps, so lambda = 1 never forgets anything.'], ...
        forgetting_lambda);
end

n_snapshots = size(snapshot_block, 2);
if n_snapshots < 1
    error('sample_covariance:EmptyBlock', 'snapshot_block contains no snapshots.');
end

% ────────────────────────── THIS STEP'S ESTIMATE ──────────────────

% Outer product averaged over the block. Rank is at most K, which is why K must
% exceed n_elements for this to be invertible on its own (see simulate_snapshots).
block_covariance = (snapshot_block * snapshot_block') / n_snapshots;

% ────────────────────────── BLEND WITH THE PAST ───────────────────

if isempty(previous_covariance)
    % First step: there is no past. Starting from zero instead would bias the
    % covariance low for the first horizon and make the early weights nonsense.
    covariance = block_covariance;
else
    covariance = forgetting_lambda * previous_covariance ...
                 + (1 - forgetting_lambda) * block_covariance;
end

% A covariance is Hermitian by construction; floating-point accumulation over
% many steps slowly breaks that, and an eigendecomposition of a matrix that is
% not quite Hermitian returns complex eigenvalues that are hard to interpret.
% Symmetrising costs nothing and keeps every downstream result real where it
% should be real.
covariance = (covariance + covariance') / 2;

% Never let a silently broken run look like a good one: a non-finite covariance
% produces a singular solve, garbage weights and a completed run with no error.
if ~all(isfinite(covariance(:)))
    error('sample_covariance:NonFinite', ...
        ['The covariance contains non-finite values. Something upstream ' ...
         'produced Inf or NaN -- check the snapshot block and the scenario ' ...
         'power levels.']);
end
end
