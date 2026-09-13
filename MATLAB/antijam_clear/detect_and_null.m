function [weights, reason] = detect_and_null(array, jammer_state)
% ══════════════════════════════════════════════════════════════════
% DETECT_AND_NULL
% Approach 1: hold the signal, put a zero on the jammer. Pure geometry.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   [weights, reason] = DETECT_AND_NULL(array, jammer_state)
%
%   Two directions, two demands, one closed-form answer:
%
%       C = [e_s , e_j]        the directions you care about
%       g = [1 ; 0]            what you want in each: keep, kill
%       w = C * inv(C' * C) * g
%
%   That is the minimum-norm solution of C' * w = g. Among all weight vectors
%   that hold unit gain on the signal AND place an exact zero on the jammer, it
%   picks the one with the smallest norm -- and since the noise that reaches the
%   output is sigma_n^2 * ||w||^2, "smallest weights that satisfy the
%   constraints" is not an aesthetic preference, it is the best of them.
%
%   WHAT IS REMARKABLE ABOUT THIS, AND WORTH SAYING OUT LOUD: there is no
%   covariance in it. No data, no inverse of an estimated quantity, no diagonal
%   loading, nothing to converge. It needs only the two directions, and the
%   whole of it can be drawn on a whiteboard. That is what makes Approach 1
%   genuinely independent of Approach 2 rather than a variation on it -- the two
%   share the detection-free solver's spine nowhere at all.
%
%   WHAT IT GIVES UP, which is measured rather than asserted. It shapes the
%   pattern at exactly two points and is indifferent to everything else:
%
%       - it ignores the noise floor beyond minimising ||w||;
%       - it ignores how strong the jammer actually is. A jammer 20 dB down
%         receives the same infinitely deep null as one 40 dB up, and the
%         degrees of freedom are spent either way;
%       - it is only as good as the estimated angle. Approach 2 never needs an
%         angle and so cannot be wrong about one.
%
%   GRACEFUL DEGRADATION -- three ways this correctly declines to null:
%
%       1. No jammer detected. Return the quiescent beam.
%       2. No degrees of freedom (a one-element array). One weight cannot
%          satisfy two constraints. Return the quiescent beam and say so.
%       3. The two directions are nearly the same VECTOR to this array. Then the
%          constraints fight each other and ||w|| explodes. This is not just the
%          jammer being angularly close -- the guard sector already prevents
%          that. On a mirror-ambiguous array a jammer at 180 - theta_s is a full
%          hundred degrees away and STILL has the same steering vector as the
%          signal. Nulling it would cancel the signal exactly. Refusing is the
%          only correct answer, and detecting it needs the vectors, not the
%          angles.
%
%   Inputs:
%       array        : struct from make_array.
%       jammer_state : struct from detect_jammer. Only its is_present,
%                      theta_deg and phi_deg fields are used -- where to null,
%                      not how the decision was reached.
%
%   Outputs:
%       weights : (n_elements x 1) complex, unit gain on the signal.
%                 Units: dimensionless.
%       reason  : text describing what was done, for reporting.

% ────────────────────────── CONSTANTS ─────────────────────────────

% Normalised inner product above which the signal and jammer directions are the
% same vector as far as this array is concerned. The 2x2 constraint matrix then
% has condition number about (1 + c)/(1 - c), so c = 0.99 already means a
% hundredfold amplification of the weights; beyond it the "null" is numerical
% noise and the signal goes with it.
CONSTRAINT_COHERENCE_LIMIT = 0.99;

profile = array.profile;

signal_steering = steering_vector(array.element_patterns, array.theta_deg, ...
                                  array.phi_deg, profile.signal_theta_deg, ...
                                  profile.signal_phi_deg);

% ────────────────────────── CASE 1: NOTHING TO NULL ───────────────

if ~jammer_state.is_present || ~isfinite(jammer_state.theta_deg)
    weights = quiescent_weights(signal_steering);

    % When the detector could not run at all -- too few elements for a noise
    % subspace -- it knows why, and that explanation is far more useful than
    % "no jammer detected". Pass it through rather than restating the symptom.
    if isfield(jammer_state, 'is_feasible') && ~jammer_state.is_feasible ...
            && isfield(jammer_state, 'reason') && ~isempty(jammer_state.reason)
        reason = sprintf('%s Quiescent beam.', jammer_state.reason);
    else
        reason = 'No jammer detected: quiescent beam, all gain on the signal.';
    end
    return
end

% ────────────────────────── CASE 2: NO FREEDOM ────────────────────

% Via detect_jammer this branch is unreachable, because an array with no
% degrees of freedom also has too few elements for MUSIC and so never reports a
% jammer -- case 1 catches it first, with a better explanation. It is kept
% because detect_and_null takes a jammer direction from whatever supplies one,
% and a solve that cannot exist should be refused by the solver rather than
% only by its usual caller.
if profile.degrees_of_freedom < 1
    weights = quiescent_weights(signal_steering);
    reason  = sprintf( ...
        ['This array has %d element and therefore %d degrees of freedom after ' ...
         'the signal constraint. One weight cannot satisfy two constraints, so ' ...
         'no null exists at any angle. Quiescent beam.'], ...
        profile.n_elements, profile.degrees_of_freedom);
    return
end

% ────────────────────────── CASE 3: DIRECTIONS TOO ALIKE ──────────

jammer_steering = steering_vector(array.element_patterns, array.theta_deg, ...
                                  array.phi_deg, jammer_state.theta_deg, ...
                                  jammer_state.phi_deg);

norm_product = norm(signal_steering) * norm(jammer_steering);
if norm_product <= 0
    weights = quiescent_weights(signal_steering);
    reason  = ['The array radiates nothing in the estimated jammer direction, ' ...
               'so there is nothing to null there. Quiescent beam.'];
    return
end

constraint_coherence = abs(signal_steering' * jammer_steering) / norm_product;

if constraint_coherence > CONSTRAINT_COHERENCE_LIMIT
    weights = quiescent_weights(signal_steering);
    reason  = sprintf( ...
        ['The jammer direction (%.1f, %.1f) deg has the same steering vector ' ...
         'as the signal to within %.4f, so a null there would cancel the ' ...
         'wanted signal. Refusing to null; quiescent beam.'], ...
        jammer_state.theta_deg, jammer_state.phi_deg, constraint_coherence);
    return
end

% ────────────────────────── THE SOLVE ─────────────────────────────

constraint_matrix   = [signal_steering, jammer_steering];        % (n_el x 2)
constraint_response = [1; 0];                                    % keep; kill

% Minimum-norm solution of constraint_matrix' * w = constraint_response.
% The 2x2 Gram matrix is well conditioned because of the check above, so the
% backslash is a two-by-two solve and nothing more.
gram_matrix = constraint_matrix' * constraint_matrix;
weights     = constraint_matrix * (gram_matrix \ constraint_response);

if ~all(isfinite(weights))
    error('detect_and_null:NonFiniteWeights', ...
        ['The constraint solve produced non-finite weights, which means the ' ...
         'two directions were closer than the coherence check allowed. ' ...
         'Coherence was %.6f.'], constraint_coherence);
end

reason = sprintf( ...
    'Null placed at (%.1f, %.1f) deg; jammer classified %s.', ...
    jammer_state.theta_deg, jammer_state.phi_deg, jammer_state.behaviour);
end
