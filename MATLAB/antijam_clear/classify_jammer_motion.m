function motion = classify_jammer_motion(theta_history_deg, phi_history_deg, ...
                                         presence_history, profile, ...
                                         covariance_horizon_steps, estimate_lag_steps)
% ══════════════════════════════════════════════════════════════════
% CLASSIFY_JAMMER_MOTION
% What is the jammer doing, and therefore where should the null go?
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   motion = CLASSIFY_JAMMER_MOTION(theta_history_deg, phi_history_deg, ...
%                                   presence_history, profile, ...
%                                   covariance_horizon_steps, estimate_lag_steps)
%
%   estimate_jammer_angle answers "where is it, right now". This answers "what
%   has it been doing", which is a different question with different evidence --
%   a history rather than a single covariance. Keeping the two apart is
%   deliberate: fusing them is how the previous implementation grew a single
%   485-line function that could not be tested or explained in pieces.
%
%   THE DECISION, IN THE ORDER IT RUNS.
%
%   0. FOLD FIRST, if the array is mirror-ambiguous. On such an array e(theta)
%      and e(180 - theta) are the same vector, so MUSIC returns one or the other
%      arbitrarily and an unfolded history looks like a jammer teleporting a
%      hundred degrees between steps. Folding every angle into one half-space
%      before looking at motion is what stops a stationary jammer being
%      classified as violently drifting. Nulling is unaffected either way --
%      a null at the mirror IS the null at the truth, because the steering
%      vectors are identical.
%
%   1. SEEN INTERMITTENTLY  -> 'onoff'.   Hold the last known position.
%   2. SEEN MOVING STEADILY -> 'drifting'. Extrapolate ahead (see below).
%   3. SEEN CONSISTENTLY    -> 'steady'.   Average the window down.
%   4. NOT SEEN AT ALL      -> 'absent'.   No null to place.
%
%   Each behaviour gets the aiming rule that suits it, and each rule is one
%   sentence: if it is not moving, average out the jitter; if it comes and goes,
%   hold where you last saw it; if it is moving, lead it.
%
%   THE LEAD, AND WHY IT IS NOT A TUNING KNOB.
%
%       aim = latest estimate + (measured drift rate) x (estimate lag)
%
%   Both quantities on the right are already known: the rate is the slope this
%   function fits, and the lag follows from lambda. Nothing is fitted to
%   performance and there is no coefficient to choose. If lambda changes, the
%   lead changes with it automatically.
%
%   WHICH "HORIZON" -- THIS DISTINCTION IS WORTH ONE dB. The covariance weights
%   a block from k steps ago by (1 - lambda) * lambda^k, and that decaying
%   weighting has two different one-number summaries:
%
%       1/(1 - lambda)       = 10 steps at lambda 0.90 -- the EFFECTIVE WINDOW
%                              LENGTH, i.e. how much smoothing this is
%                              equivalent to. Sets how long the array remembers,
%                              and so sets the motion window below.
%
%       lambda/(1 - lambda)  =  9 steps at lambda 0.90 -- the MEAN AGE of the
%                              data in that weighted average. The beamscan
%                              reports where the jammer was ON AVERAGE, and that
%                              average is this many steps old. THIS is the lag,
%                              and therefore this is what the lead must use.
%
%   They differ by exactly 1, since 1/(1-L) - L/(1-L) = 1. Using the window
%   length for the lead over-shoots by one step's worth of motion every step.
%   Measured directly: at drift rates of 0.10 and 0.20 deg/step the lag divided
%   by the rate comes out at 9.00 steps exactly. Correcting it is worth +0.4 to
%   +1.3 dB across the qualified drift range and takes the demo scenario from
%   10.8 dB / score 98 to 11.1 dB / score 99, with the mean absolute aim error
%   falling from 3.96 to 0.75 degrees.
%
%   This is NOT the constant-velocity Kalman predictor of the previous
%   implementation, which estimated angular velocity with a proper filter and
%   used a longer lead of 14 steps (10 for the covariance plus 4 for
%   application and weight-smoothing delay, neither of which exists here). That
%   predictor was worth +35 points on drifting jammers. This is the arithmetic
%   core of the same idea with none of the filtering, and its benefit is
%   measured rather than assumed -- see README.md, "Deliberate omissions".
%
%   Inputs:
%       theta_history_deg       : (n_steps x 1) estimated elevations, NaN on
%                                 steps where the jammer was not detected.
%                                 Units: degrees.
%       phi_history_deg         : (n_steps x 1) estimated azimuths, NaN as above.
%                                 Units: degrees.
%       presence_history        : (n_steps x 1) logical, jammer detected.
%       profile                 : struct from array_profile. Supplies
%                                 is_mirror_ambiguous.
%       covariance_horizon_steps: 1/(1 - forgetting_lambda), the effective
%                                 window length. Sets the motion window.
%                                 Units: steps.
%       estimate_lag_steps      : forgetting_lambda/(1 - forgetting_lambda), the
%                                 mean age of the covariance. Sets the lead.
%                                 Units: steps.
%
%   Outputs:
%       motion : struct with fields
%           behaviour            : 'steady' | 'onoff' | 'drifting' | 'absent'.
%           aim_theta_deg        : where to place the null. NaN if absent.
%           aim_phi_deg          : ditto. Units: degrees.
%           drift_rate_deg_step  : fitted elevation rate. Units: deg/step.
%           lead_deg             : the extrapolation applied. Units: degrees.
%           duty_cycle           : fraction of the window the jammer was seen.
%           n_observations       : detections inside the window.
%           reason               : one-line explanation of the classification.

% ────────────────────────── CONSTANTS ─────────────────────────────

% How far back to look, in covariance horizons. Three horizons is the shortest
% window in which the covariance has fully refreshed more than twice, so a
% change of state seen inside it is a real change and not leftover memory.
MOTION_WINDOW_HORIZONS = 3;

% Seen on fewer than this fraction of the window's steps means the jammer is
% switching rather than steady. A genuinely steady strong jammer reads 1.0;
% 0.9 allows one missed detection in ten before the verdict changes.
STEADY_DUTY_THRESHOLD = 0.9;

% Angle estimates are quantised to the pattern grid (1 deg here), so across a
% 30-step window the smallest slope distinguishable from quantisation is about
% one grid step per window, 0.033 deg/step. Three times that is required before
% motion is called drift: below it, the total movement across the whole window
% is under three grid cells and indistinguishable from rounding.
GRID_STEPS_TO_CALL_DRIFT = 3;
ANGLE_GRID_STEP_DEG      = 1.0;

% Elevation is physically bounded; an extrapolation may not leave the sphere.
THETA_MIN_DEG = 0.0;
THETA_MAX_DEG = 180.0;

theta_history_deg = theta_history_deg(:);
phi_history_deg   = phi_history_deg(:);
presence_history  = logical(presence_history(:));

motion = struct( ...
    'behaviour',           'absent', ...
    'aim_theta_deg',       NaN, ...
    'aim_phi_deg',         NaN, ...
    'drift_rate_deg_step', 0.0, ...
    'lead_deg',            0.0, ...
    'duty_cycle',          0.0, ...
    'n_observations',      0, ...
    'reason',              '');

% ────────────────────────── FOLD THE MIRROR ───────────────────────

% Only on arrays where the two half-spaces are genuinely indistinguishable.
% Doing it unconditionally would destroy real information on an array that CAN
% tell them apart, such as Monopoles (mirror coherence 0.19).
if profile.is_mirror_ambiguous
    theta_history_deg = fold_to_upper_half(theta_history_deg);
end

% ────────────────────────── THE WINDOW ────────────────────────────

window_length = MOTION_WINDOW_HORIZONS * covariance_horizon_steps;
n_steps       = numel(theta_history_deg);
first_step    = max(1, n_steps - window_length + 1);

window_index    = (first_step : n_steps).';
window_theta    = theta_history_deg(window_index);
window_phi      = phi_history_deg(window_index);
window_presence = presence_history(window_index);

motion.duty_cycle     = mean(window_presence);
motion.n_observations = sum(window_presence);

if motion.n_observations == 0
    motion.reason = sprintf( ...
        'No detection in the last %d steps: no jammer to null.', numel(window_index));
    return
end

% Only the steps where something was actually seen carry angle information.
seen_step_offsets = find(window_presence);
seen_theta        = window_theta(seen_step_offsets);
seen_phi          = window_phi(seen_step_offsets);

% ────────────────────────── INTERMITTENT? ─────────────────────────

if motion.duty_cycle < STEADY_DUTY_THRESHOLD
    % Hold the last place it was actually seen. The null costs almost nothing
    % while the jammer is off, and it is already in place when it returns --
    % which is why a stale null is the right answer here rather than a
    % second-best one.
    motion.behaviour     = 'onoff';
    motion.aim_theta_deg = seen_theta(end);
    motion.aim_phi_deg   = seen_phi(end);
    motion.reason        = sprintf( ...
        ['Detected on %.0f%% of the last %d steps: the jammer is switching. ' ...
         'Holding the null at the last known position.'], ...
        100 * motion.duty_cycle, numel(window_index));
    return
end

% ────────────────────────── MOVING OR STILL? ──────────────────────

% Straight-line fit of elevation against step number over the detections in the
% window. Slope is the angular rate.
if numel(seen_theta) >= 2
    line_coefficients         = polyfit(seen_step_offsets, seen_theta, 1);
    motion.drift_rate_deg_step = line_coefficients(1);
else
    motion.drift_rate_deg_step = 0.0;
end

drift_threshold_deg_step = GRID_STEPS_TO_CALL_DRIFT * ANGLE_GRID_STEP_DEG ...
                           / max(numel(window_index), 1);

if abs(motion.drift_rate_deg_step) > drift_threshold_deg_step
    % Moving. Aim ahead by exactly the covariance's own lag -- see the header.
    % The lead uses the MEAN AGE of the covariance, not its window length --
    % see the header. These differ by exactly one step, and using the wrong one
    % over-leads by one step's worth of motion on every step.
    motion.behaviour = 'drifting';
    motion.lead_deg  = motion.drift_rate_deg_step * estimate_lag_steps;

    % CAP THE LEAD AT THE GUARD SECTOR -- half a beamwidth. The lead multiplies
    % the fitted rate by the mean age, so it amplifies any error in that rate
    % ninefold. That is a good trade while the rate is well measured, and a bad
    % one once it is not: past roughly one beamwidth of travel per horizon the
    % beamscan peak itself smears, the fitted slope gets noisy, and the
    % extrapolation throws the null further off than no extrapolation at all.
    % Beyond half a beamwidth ahead there is no measurement bearing on where
    % the jammer will be, so the cap says: do not aim further than you can see.
    % It reuses guard_deg and introduces no new constant.
    maximum_lead_deg = profile.guard_deg;
    if abs(motion.lead_deg) > maximum_lead_deg
        motion.lead_deg = sign(motion.lead_deg) * maximum_lead_deg;
    end

    aim_theta_deg        = seen_theta(end) + motion.lead_deg;
    motion.aim_theta_deg = min(max(aim_theta_deg, THETA_MIN_DEG), THETA_MAX_DEG);
    motion.aim_phi_deg   = seen_phi(end);

    motion.reason = sprintf( ...
        ['Moving at %.2f deg/step. Aiming %.1f deg ahead of the covariance ' ...
         'estimate, which describes where the jammer was %d steps ago.'], ...
        motion.drift_rate_deg_step, motion.lead_deg, estimate_lag_steps);
else
    % Still. Averaging the window removes the grid quantisation jitter, and
    % the median is used so a single bad estimate cannot drag the null off.
    motion.behaviour     = 'steady';
    motion.aim_theta_deg = median(seen_theta);
    motion.aim_phi_deg   = median(seen_phi);

    motion.reason = sprintf( ...
        ['Stationary: fitted rate %.2f deg/step from %d observation(s), below ' ...
         'the %.2f deg/step needed to call drift. Nulling the window median.'], ...
        motion.drift_rate_deg_step, numel(seen_theta), drift_threshold_deg_step);
end
end


% ────────────────────────── HELPERS ───────────────────────────────

function folded_theta_deg = fold_to_upper_half(theta_deg)
% Map elevation into [0, 90] by reflecting about 90 degrees.
%
% Only valid where e(theta) and e(180 - theta) are the same steering vector, so
% that the two angles are physically indistinguishable to the array and either
% may be used to place the null. The information discarded -- which side the
% jammer is on -- is information the array never had.
folded_theta_deg = theta_deg;
is_lower_half    = theta_deg > 90.0;
folded_theta_deg(is_lower_half) = 180.0 - theta_deg(is_lower_half);
end
