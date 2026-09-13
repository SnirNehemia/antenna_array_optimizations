function weights = max_sinr_weights(covariance, array, loading_factor)
% ══════════════════════════════════════════════════════════════════
% MAX_SINR_WEIGHTS
% Approach 2: maximise SINR directly, never asking where the jammer is.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   weights = MAX_SINR_WEIGHTS(covariance, array, loading_factor)
%
%   One line of algebra:
%
%       R_loaded = R + delta * I
%       w        = R_loaded \ e_s  /  (e_s' * (R_loaded \ e_s))
%
%   WHY THIS MAXIMISES SINR. It is the minimum-variance solution: among all
%   weight vectors with unit gain on the signal, it is the one that minimises
%   the total output power w' R w. The signal's contribution is pinned at unit
%   gain by the constraint, so minimising the total minimises EVERYTHING ELSE --
%   jammer and noise together, in whatever proportion they happen to arrive.
%
%   It never asks where the jammer is. It does not know, and cannot be told. It
%   simply refuses to pass power that is not coming from the signal direction,
%   and a jammer is by definition power that is not coming from there. That is
%   the whole of Approach 2, and it is why it cannot be wrong about an angle:
%   it never estimates one.
%
%   WHY DIAGONAL LOADING IS STRUCTURAL AND NOT A KNOB. The wanted signal is
%   inside R -- a real receiver cannot remove what it is trying to receive.
%   That makes this MPDR rather than MVDR, and it has a sharp consequence: if
%   the assumed signal direction e_s is even slightly wrong, the solver sees
%   wanted signal arriving from a direction it is not protecting, correctly
%   identifies it as "power that is not the signal", and cancels YOUR OWN
%   TRANSMISSION. The stronger your signal, the more of it there is to cancel,
%   so this failure gets worse exactly when things should be getting easier.
%
%   Adding delta * I raises the noise floor the solver believes in, which caps
%   how much of any single direction it is willing to cancel -- including
%   yours. Without it, a steering error of one grid cell is enough to do real
%   damage. This is not a performance tweak; it is what makes the method work
%   at all on measured patterns, where the true steering vector is never
%   exactly the assumed one.
%
%   HOW DELTA IS SET.
%
%       delta = loading_factor * (estimated noise floor)
%
%   Scaled to the NOISE floor, not to trace(R). The trace is dominated by the
%   jammer, so trace scaling would raise your own loading whenever the jammer
%   got stronger -- blunting the very null you need against it. The noise floor
%   is the one part of the covariance that does not move when the jammer does.
%   See noise_floor_power.
%
%   loading_factor is the single genuine tuning parameter in this folder. It is
%   chosen by measuring performance against STEERING MISMATCH, not against null
%   depth: null depth always prefers less loading, and always being wrong in
%   that direction is how self-cancellation gets shipped.
%
%   Inputs:
%       covariance     : (n_elements x n_elements) Hermitian, from
%                        sample_covariance. Units: power.
%       array          : struct from make_array. Supplies the patterns, grids
%                        and the assumed signal direction.
%       loading_factor : delta as a multiple of the estimated noise floor.
%                        Units: dimensionless. Must be > 0.
%
%   Outputs:
%       weights : (n_elements x 1) complex, unit gain on the assumed signal
%                 direction. Units: dimensionless.

if ~(loading_factor > 0)
    error('max_sinr_weights:BadLoading', ...
        ['loading_factor must be positive; got %g. Zero loading is not an ' ...
         'option here: the wanted signal is inside the covariance, so an ' ...
         'unloaded solve cancels it whenever the steering vector is imperfect.'], ...
        loading_factor);
end

profile = array.profile;

signal_steering = steering_vector(array.element_patterns, array.theta_deg, ...
                                  array.phi_deg, profile.signal_theta_deg, ...
                                  profile.signal_phi_deg);

% ────────────────────────── LOAD THE DIAGONAL ─────────────────────

noise_power = noise_floor_power(covariance);
delta       = loading_factor * noise_power;

n_elements       = size(covariance, 1);
loaded_covariance = covariance + delta * eye(n_elements);

% ────────────────────────── THE SOLVE ─────────────────────────────

% Solved rather than inverted: forming inv(R) explicitly is both slower and
% less accurate, and nothing here needs the inverse itself.
whitened_steering = loaded_covariance \ signal_steering;

% Normalise for unit gain on the signal: w' * e_s = 1.
normalisation = signal_steering' * whitened_steering;
if abs(normalisation) <= 0
    error('max_sinr_weights:DegenerateSolve', ...
        ['The normalisation e_s'' R^-1 e_s came out zero, which means the ' ...
         'array has no response in the assumed signal direction. Check the ' ...
         'target angle against the element patterns.']);
end

% The conjugate keeps the convention y = w'*x, under which w'*e_s is exactly 1.
weights = whitened_steering / conj(normalisation);

% Never let a silently broken solve look like a good run.
if ~all(isfinite(weights))
    error('max_sinr_weights:NonFiniteWeights', ...
        ['The loaded solve produced non-finite weights. The covariance is ' ...
         'singular despite loading of %g -- check that it is finite and that ' ...
         'the noise floor estimate is sane.'], delta);
end
end
