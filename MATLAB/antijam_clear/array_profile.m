function profile = array_profile(element_patterns, theta_deg, phi_deg, ...
                                 signal_theta_deg, signal_phi_deg)
% ══════════════════════════════════════════════════════════════════
% ARRAY_PROFILE
% Measure the array once, before adapting, and derive the run's constants.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   profile = ARRAY_PROFILE(element_patterns, theta_deg, phi_deg, ...
%                           signal_theta_deg, signal_phi_deg)
%
%   Three facts that the rest of the run depends on, each measured from the
%   array itself rather than assumed:
%
%   1. GUARD SECTOR (profile.guard_deg). A jammer angularly closer to the target
%      than the beamwidth is INSIDE the main lobe. Nulling it cancels the wanted
%      signal along with it, so such a geometry is outside what this method
%      claims to do, and it must be recognised rather than scored as a failure.
%      The guard is derived per (array, target) as half the wider 3 dB
%      beamwidth. It is NOT a global constant: measured across the arrays here
%      the true value spans roughly 16 to 72 degrees, so any single number is
%      wrong on nearly every array.
%
%   2. MIRROR AMBIGUITY (profile.mirror_coherence). Some arrays -- planar ones
%      especially -- have e(theta, phi) and e(180 - theta, phi) essentially
%      identical, to parts in 1e6. This is BENIGN for nulling: the two steering
%      vectors are the same vector, so a null placed at the mirror IS the null
%      at the truth. It is FATAL for tracking: an angle estimate that flips
%      between theta and 180 - theta between steps looks like violent motion,
%      and a stationary jammer gets classified as a drifting one. Flagged here;
%      folded in classify_jammer_motion.
%
%      Note the pair tested is the MIRROR (180 - theta, phi), not the antipode
%      (180 - theta, phi + 180). Confusing the two produces a confident and
%      wrong conclusion that the ambiguity does not exist.
%
%   3. DEGREES OF FREEDOM (profile.degrees_of_freedom). An n-element array has n
%      complex weights. One is spent holding gain on the signal; the remainder
%      are what is available to place nulls with. At n_elements = 1 that count
%      is zero and the array provably cannot null anything -- the correct output
%      is then the quiescent beam and a stated reason, not an exception. This is
%      how the method demonstrates that it knows its own limits.
%
%   4. TARGET VISIBILITY (profile.target_visibility_db). How much weaker the
%      array's best possible beam is at the TARGET than in the array's best
%      direction anywhere on the sphere. If the target sits in a deep pattern
%      null the array cannot receive its own signal there, so no beamformer --
%      this one, the previous one, or the perfect-knowledge oracle -- can
%      achieve anything, and a score computed for that geometry is meaningless
%      rather than bad. Such a pair is NOT A TEST and must be reported as such.
%
%      Measured on the arrays here this quantity is strongly bimodal: genuine
%      geometries sit 7.6 to 15.7 dB below peak, while patch_back2back at
%      (30, 0) sits 113 dB below -- back-to-back patches radiate no co-pol in
%      that direction at all. The threshold is placed in the middle of that
%      ~97 dB gap, so its exact value cannot change any conclusion.
%
%   Inputs:
%       element_patterns : (n_elements x n_theta x n_phi) complex. Units: V/m.
%       theta_deg        : (n_theta x 1) elevation grid. Units: degrees.
%       phi_deg          : (n_phi x 1)   azimuth grid.   Units: degrees.
%       signal_theta_deg : wanted-signal elevation. Units: degrees.
%       signal_phi_deg   : wanted-signal azimuth.   Units: degrees.
%
%   Outputs:
%       profile : struct with fields
%           signal_theta_deg      : the target this profile was measured for.
%           signal_phi_deg        : ditto. Every field below is a property of
%                                   the (array, target) PAIR, not of the array
%                                   alone, so the target travels with them.
%           n_elements            : element count.
%           degrees_of_freedom    : n_elements - 1, nulls placeable.
%           can_null              : logical, degrees_of_freedom >= 1.
%           quiescent_dbi         : directivity of the plain beam at the target.
%                                   Units: dBi.
%           hpbw_theta_deg        : 3 dB beamwidth in the theta cut. Units: deg.
%           hpbw_phi_deg          : 3 dB beamwidth in the phi cut.   Units: deg.
%           guard_deg             : derived guard half-sector. Units: degrees.
%           mirror_coherence      : 0..1 similarity of e(theta) and e(180-theta).
%           is_mirror_ambiguous   : logical, mirror_coherence > 0.99.
%           target_visibility_db  : target strength relative to the array's best
%                                   direction, <= 0. Units: dB.
%           is_target_illuminated : logical, target_visibility_db > -40.

% ────────────────────────── CONSTANTS ─────────────────────────────

% Above this normalised inner product the two mirror directions are the same
% steering vector to within numerical noise, and any direction TRACKED over
% time must be folded into one half-space. 0.99 is far above every genuinely
% distinct pair measured here and far below the ambiguous ones.
MIRROR_COHERENCE_THRESHOLD = 0.99;

% The guard half-sector as a fraction of the wider beamwidth. 0.5 says: closer
% to the target than half a beamwidth counts as inside the main lobe. This is a
% stated judgement, not a derivation -- it lives here so it can be defended or
% changed in exactly one place.
GUARD_FRACTION_OF_BEAMWIDTH = 0.5;

% A target more than this far below the array's best direction is in a pattern
% null: the array cannot receive its own signal there, so the geometry is not a
% test. See note 4 above -- the measured values are bimodal with a ~97 dB gap,
% and this sits inside the gap, so the precise number is not load-bearing.
TARGET_VISIBILITY_FLOOR_DB = -40.0;

n_elements = size(element_patterns, 1);
theta_deg  = theta_deg(:);
phi_deg    = phi_deg(:);

% ────────────────────────── THE QUIESCENT BEAM ────────────────────

% Beamwidth and directivity are properties of the beam you actually form when
% you are not adapting, so both are measured on the quiescent (matched-filter)
% beam steered at the target.
signal_steering = steering_vector(element_patterns, theta_deg, phi_deg, ...
                                  signal_theta_deg, signal_phi_deg);
weights = quiescent_weights(signal_steering);

% compute_array_factor forms sum_n w_n * E_n, which is w.' * e, whereas this
% folder's convention is y = w' * x. conj() reconciles them -- see
% steering_vector for the full note.
array_factor         = compute_array_factor(conj(weights), element_patterns);
directivity_dbi_grid = compute_directivity_dbi_grid(array_factor, theta_deg, phi_deg, []);

theta_index = nearest_index(theta_deg, signal_theta_deg);
phi_index   = nearest_index(phi_deg,   signal_phi_deg);

quiescent_dbi = directivity_dbi_grid(theta_index, phi_index);

% ────────────────────────── BEAMWIDTHS ────────────────────────────

theta_step_deg = grid_step(theta_deg);
phi_step_deg   = grid_step(phi_deg);

% Theta cut at the signal azimuth. The cut at phi + 180 completes the great
% circle, so compute_hpbw can follow a beam that runs over the pole.
opposite_phi_index = nearest_index(phi_deg, wrap_phi(signal_phi_deg + 180.0, phi_deg));
theta_cut_dbi      = directivity_dbi_grid(:, phi_index);
opposite_cut_dbi   = directivity_dbi_grid(:, opposite_phi_index);

% compute_hpbw takes a 0-BASED peak index (it mirrors the Python original).
theta_full_span_deg = (numel(theta_deg) - 1) * theta_step_deg;
hpbw_theta_deg = full_span_if_no_crossing( ...
    compute_hpbw(theta_cut_dbi, theta_index - 1, opposite_cut_dbi) * theta_step_deg, ...
    theta_full_span_deg);

% Phi cut at the signal elevation. Phi is circular and compute_hpbw scans
% linearly, so the cut is rotated to put the peak in the middle first --
% otherwise a beam straddling phi = 0 reads as two half-beams.
phi_cut_dbi     = directivity_dbi_grid(theta_index, :).';
shift_amount    = floor(numel(phi_deg) / 2) - (phi_index - 1);
phi_cut_centred = circshift(phi_cut_dbi, shift_amount);
centred_peak_index_0based = floor(numel(phi_deg) / 2);

phi_full_span_deg = numel(phi_deg) * phi_step_deg;
hpbw_phi_deg = full_span_if_no_crossing( ...
    compute_hpbw(phi_cut_centred, centred_peak_index_0based, []) * phi_step_deg, ...
    phi_full_span_deg);

guard_deg = GUARD_FRACTION_OF_BEAMWIDTH * max(hpbw_theta_deg, hpbw_phi_deg);

% ────────────────────────── MIRROR AMBIGUITY ──────────────────────

mirror_steering = steering_vector(element_patterns, theta_deg, phi_deg, ...
                                  180.0 - signal_theta_deg, signal_phi_deg);

% Normalised inner product: 1 means the two directions are indistinguishable to
% the array, 0 means they are orthogonal.
norm_product = norm(signal_steering) * norm(mirror_steering);
if norm_product <= 0
    mirror_coherence = 0.0;
else
    mirror_coherence = abs(signal_steering' * mirror_steering) / norm_product;
end

% ────────────────────────── TARGET VISIBILITY ─────────────────────

% The matched filter's output power in a direction is ||e||^2 = sum_n |E_n|^2,
% so the best the array can ever do in ANY direction is the largest value of
% that sum over the grid. The target's share of it is the visibility.
power_over_grid    = squeeze(sum(abs(element_patterns) .^ 2, 1));   % (n_theta x n_phi)
best_power         = max(power_over_grid(:));
target_power       = real(signal_steering' * signal_steering);

if best_power <= 0 || target_power <= 0
    target_visibility_db = -Inf;
else
    target_visibility_db = 10 * log10(target_power / best_power);
end

% ────────────────────────── ASSEMBLE ──────────────────────────────

profile = struct( ...
    'signal_theta_deg',    signal_theta_deg, ...
    'signal_phi_deg',      signal_phi_deg, ...
    'n_elements',          n_elements, ...
    'degrees_of_freedom',  n_elements - 1, ...
    'can_null',            n_elements - 1 >= 1, ...
    'quiescent_dbi',       quiescent_dbi, ...
    'hpbw_theta_deg',      hpbw_theta_deg, ...
    'hpbw_phi_deg',        hpbw_phi_deg, ...
    'guard_deg',           guard_deg, ...
    'mirror_coherence',    mirror_coherence, ...
    'is_mirror_ambiguous', mirror_coherence > MIRROR_COHERENCE_THRESHOLD, ...
    'target_visibility_db', target_visibility_db, ...
    'is_target_illuminated', target_visibility_db > TARGET_VISIBILITY_FLOOR_DB);
end


% ────────────────────────── HELPERS ───────────────────────────────

function step_deg = grid_step(grid_deg)
% Angular spacing of a uniform grid, in degrees.
if numel(grid_deg) < 2
    step_deg = 1.0;
else
    step_deg = mean(diff(grid_deg));
end
end


function phi_deg_wrapped = wrap_phi(phi_value_deg, phi_deg)
% Wrap an azimuth into the span covered by the phi grid.
span_deg = numel(phi_deg) * grid_step(phi_deg);
phi_deg_wrapped = mod(phi_value_deg - phi_deg(1), span_deg) + phi_deg(1);
end


function width_deg = full_span_if_no_crossing(width_deg, full_span_deg)
% compute_hpbw returns 0 when it never meets the -3 dB level, which means the
% beam never falls 3 dB below its peak anywhere in the cut. That is a real
% answer -- the pattern is effectively omnidirectional in this cut -- and the
% honest beamwidth is then the whole span, not zero.
if ~isfinite(width_deg) || width_deg <= 0
    width_deg = full_span_deg;
end
end
