function scenario = make_scenario(profile, behaviour, signal_theta_deg, signal_phi_deg, ...
                                  jammer_separation_deg, signal_to_noise_db, jammer_to_noise_db)
% ══════════════════════════════════════════════════════════════════
% MAKE_SCENARIO
% The ground truth: what the jammer really does, over time.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   scenario = MAKE_SCENARIO(profile, behaviour, signal_theta_deg, ...
%                            signal_phi_deg, jammer_separation_deg, ...
%                            signal_to_noise_db, jammer_to_noise_db)
%
%   THE ONE RULE OF THIS FILE. This struct is the only place the true jammer
%   direction exists. It is handed to exactly two things:
%
%       simulate_snapshots  -- which turns truth into data, and
%       oracle_weights / output_sinr_db -- which score against truth.
%
%   It is NEVER passed to detect_jammer, detect_and_null or max_sinr_weights.
%   An algorithm that sees the true angle will produce excellent numbers that
%   mean nothing, and the leak is invisible in the results -- they simply look
%   good. Keeping truth inside one named struct is what makes that auditable:
%   the check is a grep for 'scenario' in the algorithm files, and it must come
%   back empty.
%
%   GEOMETRY (a deliberate simplification, stated rather than hidden). The
%   jammer sits at the SAME AZIMUTH as the signal, offset in elevation by
%   jammer_separation_deg, and when it drifts it drifts in elevation, away from
%   the signal. One angular axis is enough to show every effect that matters and
%   it keeps every plot a single readable cut. A jammer offset in azimuth is a
%   change to two lines here and nothing else.
%
%   BEHAVIOURS
%       'steady' : jammer fixed in angle, transmitting throughout.
%       'onoff'  : jammer fixed in angle, switching on and off periodically.
%       'drift'  : jammer transmitting throughout, moving at a constant rate.
%
%   RUN LENGTH IS DERIVED, NOT CHOSEN. Each behaviour runs for as long as its
%   slowest estimator needs to have an opinion at all, and no fewer steps. This
%   is not fussiness: in the previous work an on/off campaign ran 4 cycles when
%   period estimation needs at least 3, so most of every run was spent before a
%   period existed. The resulting regression was read as a property of the
%   algorithm and shipped as a recommendation NOT to enable a repair that in
%   fact helps. Re-run at 8 cycles, the recommendation reversed. Run length is
%   therefore computed from the estimator's requirement and printed.
%
%   Inputs:
%       profile               : struct from array_profile. Supplies guard_deg,
%                               which decides whether this geometry is in scope.
%       behaviour             : 'steady' | 'onoff' | 'drift'.
%       signal_theta_deg      : wanted-signal elevation. Units: degrees.
%       signal_phi_deg        : wanted-signal azimuth.   Units: degrees.
%       jammer_separation_deg : jammer elevation offset from the signal, > 0.
%                               Must exceed profile.guard_deg. Units: degrees.
%       signal_to_noise_db    : wanted-signal power over noise power. Units: dB.
%       jammer_to_noise_db    : jammer power over noise power.        Units: dB.
%
%   Outputs:
%       scenario : struct with fields
%           behaviour            : the behaviour string.
%           n_steps              : number of adaptation steps (derived).
%           step_seconds         : duration of one step. Units: seconds.
%           signal_theta_deg     : scalar. Units: degrees.
%           signal_phi_deg       : scalar. Units: degrees.
%           jammer_theta_deg     : (n_steps x 1) true jammer elevation. Units: deg.
%           jammer_phi_deg       : (n_steps x 1) true jammer azimuth.   Units: deg.
%           jammer_on            : (n_steps x 1) logical, jammer transmitting.
%           onoff_period_steps   : toggle period; NaN unless behaviour is onoff.
%           drift_rate_deg_step  : drift rate; 0 unless behaviour is drift.
%           signal_power         : sigma_s^2. Units: dimensionless (noise = 1).
%           jammer_power         : sigma_j^2. Units: dimensionless.
%           noise_power          : sigma_n^2, fixed at 1 (the reference).
%           random_stream        : seeded RandStream used by simulate_snapshots.
%           description          : one-line human-readable summary.

% ────────────────────────── CONSTANTS ─────────────────────────────

% Noise power is the reference all other powers are quoted against. Fixing it
% at 1 means signal_to_noise_db and jammer_to_noise_db ARE the powers in dB.
NOISE_POWER = 1.0;

% One adaptation step per second, so a drift rate in degrees per step is also
% degrees per second and can be compared against a real target's angular rate.
STEP_SECONDS = 1.0;

% The covariance remembers 1/(1 - lambda) = 10 steps (see sample_covariance).
% Every run length below is quoted as a multiple of that horizon, so the numbers
% move together if the horizon ever changes.
COVARIANCE_HORIZON_STEPS = 10;

% The MEAN AGE of the covariance, lambda/(1 - lambda) = 9 steps. This is what a
% moving jammer's estimate actually lags by -- one step less than the window
% length above. Used only for the printed description; classify_jammer_motion
% is where it does real work.
ESTIMATE_LAG_STEPS = 9;

% A steady run needs the covariance to converge and then enough settled steps to
% average over. Six horizons gives one horizon of transient and five of measurement.
STEADY_HORIZONS = 6;

% On/off period, in covariance horizons. At exactly two horizons the jammer's
% on-time and the beamformer's memory are the same length, which is the regime
% where the reactive lag is plainly visible rather than either negligible or total.
ONOFF_PERIOD_HORIZONS = 2;

% Period estimation needs at least 3 complete cycles to have an opinion. Eight
% leaves the great majority of the run in the regime being measured rather than
% in the estimator's warm-up. This is the constant whose previous value of 4
% produced a wrong published recommendation.
ONOFF_CYCLES = 8;

% Drift rate, chosen so that the covariance's own lag is a visible fraction of a
% beamwidth rather than negligible or catastrophic: the estimate trails the
% jammer by rate * mean age = 0.5 * 9 = 4.5 degrees, about a fifth of a typical
% 25 degree beamwidth here. That lag is the effect the drift case exists to show.
DRIFT_RATE_DEG_PER_STEP = 0.5;

% Total angular travel for a drift run: far enough that the jammer leaves the
% region it started in entirely.
DRIFT_TRAVEL_DEG = 60.0;

% Fixed seed: every run in this folder is reproducible by construction.
RANDOM_SEED = 20260912;

% Elevation is physically bounded; a trajectory leaving it is a scenario bug.
THETA_MIN_DEG = 0.0;
THETA_MAX_DEG = 180.0;

% ────────────────────────── VALIDATE THE GEOMETRY ─────────────────

if jammer_separation_deg <= 0
    error('make_scenario:BadSeparation', ...
        'jammer_separation_deg must be positive; got %.1f.', jammer_separation_deg);
end

% A jammer closer to the target than the guard sector is INSIDE the main lobe.
% Nulling it cancels the wanted signal, so the geometry is outside what this
% method claims to do. Refuse it here rather than score it as a failure later.
if jammer_separation_deg <= profile.guard_deg
    error('make_scenario:JammerInsideGuard', ...
        ['Jammer separation %.1f deg is inside this array''s guard sector of ' ...
         '%.1f deg (half the wider 3 dB beamwidth, %.1f deg in theta and ' ...
         '%.1f deg in phi). A jammer that close is inside the main lobe and ' ...
         'cannot be nulled without cancelling the wanted signal. Separate them ' ...
         'further, or accept that this geometry is out of scope.'], ...
        jammer_separation_deg, profile.guard_deg, ...
        profile.hpbw_theta_deg, profile.hpbw_phi_deg);
end

if ~profile.is_target_illuminated
    error('make_scenario:TargetInPatternNull', ...
        ['The array is %.1f dB down at the target direction, which is a ' ...
         'pattern null, not a beam. It cannot receive its own signal there, ' ...
         'so no beamformer can achieve anything and any score would be ' ...
         'meaningless. This (array, target) pair is not a test.'], ...
        -profile.target_visibility_db);
end

% ────────────────────────── BUILD THE TRAJECTORY ──────────────────

onoff_period_steps  = NaN;
drift_rate_deg_step = 0.0;

switch lower(behaviour)

    case 'steady'
        n_steps          = STEADY_HORIZONS * COVARIANCE_HORIZON_STEPS;
        jammer_theta_deg = repmat(signal_theta_deg + jammer_separation_deg, n_steps, 1);
        jammer_on        = true(n_steps, 1);
        description      = sprintf('steady jammer at %.1f deg separation', ...
                                   jammer_separation_deg);

    case 'onoff'
        onoff_period_steps = ONOFF_PERIOD_HORIZONS * COVARIANCE_HORIZON_STEPS;
        n_steps            = ONOFF_CYCLES * onoff_period_steps;
        jammer_theta_deg   = repmat(signal_theta_deg + jammer_separation_deg, n_steps, 1);

        % Square wave at 50% duty: on for the first half of each period.
        step_index = (0:n_steps - 1).';
        jammer_on  = mod(step_index, onoff_period_steps) < (onoff_period_steps / 2);

        description = sprintf(['on/off jammer at %.1f deg separation, ' ...
                               'period %d steps, %d cycles'], ...
                              jammer_separation_deg, onoff_period_steps, ONOFF_CYCLES);

    case 'drift'
        drift_rate_deg_step = DRIFT_RATE_DEG_PER_STEP;
        n_steps             = round(DRIFT_TRAVEL_DEG / DRIFT_RATE_DEG_PER_STEP);

        % Moves away from the signal, so the trajectory cannot wander into the
        % guard sector partway through and silently change what is being measured.
        step_index       = (0:n_steps - 1).';
        jammer_theta_deg = signal_theta_deg + jammer_separation_deg ...
                           + drift_rate_deg_step * step_index;
        jammer_on        = true(n_steps, 1);

        description = sprintf(['drifting jammer from %.1f deg separation, ' ...
                               '%.2f deg/step (null lags by ~%.1f deg)'], ...
                              jammer_separation_deg, drift_rate_deg_step, ...
                              drift_rate_deg_step * ESTIMATE_LAG_STEPS);

    otherwise
        error('make_scenario:UnknownBehaviour', ...
            ['Unknown behaviour ''%s''. Expected ''steady'', ''onoff'' or ' ...
             '''drift''.'], behaviour);
end

jammer_phi_deg = repmat(signal_phi_deg, n_steps, 1);

% ────────────────────────── VALIDATE THE TRAJECTORY ───────────────

if any(jammer_theta_deg < THETA_MIN_DEG) || any(jammer_theta_deg > THETA_MAX_DEG)
    error('make_scenario:TrajectoryLeavesSphere', ...
        ['The jammer trajectory runs from %.1f to %.1f deg elevation, outside ' ...
         'the physical range [%.0f, %.0f]. Reduce the separation, the drift ' ...
         'rate or the run length.'], ...
        min(jammer_theta_deg), max(jammer_theta_deg), THETA_MIN_DEG, THETA_MAX_DEG);
end

separation_over_run_deg = abs(jammer_theta_deg - signal_theta_deg);
if any(separation_over_run_deg <= profile.guard_deg)
    error('make_scenario:TrajectoryEntersGuard', ...
        ['The jammer trajectory closes to %.1f deg of the signal, inside the ' ...
         '%.1f deg guard sector. Part of this run would be measuring a jammer ' ...
         'inside the main lobe, which is out of scope.'], ...
        min(separation_over_run_deg), profile.guard_deg);
end

% ────────────────────────── ASSEMBLE ──────────────────────────────

scenario = struct( ...
    'behaviour',           lower(behaviour), ...
    'n_steps',             n_steps, ...
    'step_seconds',        STEP_SECONDS, ...
    'signal_theta_deg',    signal_theta_deg, ...
    'signal_phi_deg',      signal_phi_deg, ...
    'jammer_theta_deg',    jammer_theta_deg, ...
    'jammer_phi_deg',      jammer_phi_deg, ...
    'jammer_on',           jammer_on, ...
    'onoff_period_steps',  onoff_period_steps, ...
    'drift_rate_deg_step', drift_rate_deg_step, ...
    'signal_power',        NOISE_POWER * 10 ^ (signal_to_noise_db / 10), ...
    'jammer_power',        NOISE_POWER * 10 ^ (jammer_to_noise_db / 10), ...
    'noise_power',         NOISE_POWER, ...
    'random_stream',       RandStream('mt19937ar', 'Seed', RANDOM_SEED), ...
    'description',         description);
end
