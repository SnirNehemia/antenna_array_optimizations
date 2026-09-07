function [pred, state] = adapt_cv_update(state, theta_meas_deg, phi_meas_deg, has_meas)
% ADAPT_CV_UPDATE  One constant-velocity Kalman step on the jammer's (theta, phi).
%
%   [pred, state] = ADAPT_CV_UPDATE(state, theta_meas_deg, phi_meas_deg, has_meas)
%
%   Advances the CV track by one closed-loop step and returns the angle the
%   null should be aimed at. Called once per step by adapt_predict_update with
%   the MUSIC estimate for this step; has_meas = false on steps where MUSIC
%   declared the jammer absent, in which case the filter coasts on the motion
%   model (which is the entire point of having one).
%
%   Sequence per step: predict -> fold the measurement to the nearer mirror
%   branch if the array is degenerate -> gate -> update -> extrapolate
%   cv_lead_steps ahead for the caller.
%
%   THE LEAD is what actually buys the performance. The applied null is late
%   by the covariance horizon only for the reactive tracker; for this path the
%   remaining lag is the one-step application delay plus the
%   weight_smoothing_mu first-order lag. lead_steps is tuned against that, not
%   derived — see the sweep in docs/notes.md [P12b].
%
%   Inputs:
%       state          : struct from adapt_cv_init (advanced state returned).
%       theta_meas_deg : MUSIC elevation estimate this step [deg]. Ignored
%                        when has_meas is false.
%       phi_meas_deg   : MUSIC azimuth estimate this step [deg]. Ignored when
%                        has_meas is false.
%       has_meas       : logical. False = jammer not detected this step; the
%                        filter coasts.
%
%   Outputs:
%       pred  : struct with fields
%           valid          : logical. True only when the track has been
%                            accepted for min_track_steps AND the estimated
%                            angular speed exceeds min_speed_deg_s. The caller
%                            must keep its existing branch when this is false
%                            — that is what keeps STATIC and on/off runs
%                            byte-identical to the pre-CV behaviour.
%           theta_deg      : predicted elevation lead_steps ahead [deg].
%           phi_deg        : predicted azimuth lead_steps ahead [deg].
%           speed_deg_s    : estimated angular speed [deg/s].
%           theta_now_deg  : filtered elevation at the CURRENT step [deg]
%                            (diagnostics: separates filtering from leading).
%           phi_now_deg    : filtered azimuth at the current step [deg].
%           accepted       : logical, this step's measurement passed the gate.
%       state : updated filter state.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P12b].

cfg = state.cfg;
pred = struct('valid', false, 'theta_deg', NaN, 'phi_deg', NaN, ...
    'speed_deg_s', NaN, 'theta_now_deg', NaN, 'phi_now_deg', NaN, ...
    'accepted', false);

% ────────────────────────── INITIALIZE ON FIRST FIX ───────────────
if isempty(state.x)
    if ~has_meas || ~isfinite(theta_meas_deg) || ~isfinite(phi_meas_deg)
        return                                     % nothing to start from yet
    end
    % Zero initial velocity, with a covariance wide enough in the velocity
    % states that the first few measurements — not this guess — determine it.
    state.x = [theta_meas_deg; 0; wrap360(phi_meas_deg); 0];
    state.P = diag([cfg.r_meas_deg^2, (10 * cfg.q_accel_deg_s2 + 1)^2, ...
                    cfg.r_meas_deg^2, (10 * cfg.q_accel_deg_s2 + 1)^2]);
    state.n_track = 1;
    state.n_miss  = 0;
    pred.accepted = true;
    return
end

% ────────────────────────── PREDICT ───────────────────────────────
x_pred = state.F * state.x;
P_pred = state.F * state.P * state.F' + state.Q;
x_pred(3) = wrap360(x_pred(3));

% ────────────────────────── MEASUREMENT UPDATE ────────────────────
accepted = false;
if has_meas && isfinite(theta_meas_deg) && isfinite(phi_meas_deg)
    [z, ~] = fold_measurement(state, x_pred, theta_meas_deg, phi_meas_deg);

    % Innovation, with the azimuth residual wrapped to [-180, 180). Without
    % the wrap a jammer crossing phi = 0 produces a 359 deg residual and the
    % gate throws away a perfectly good measurement.
    y = [z(1) - x_pred(1); wrap180(z(2) - x_pred(3))];

    if hypot(y(1), y(2)) <= cfg.gate_deg
        S = state.H * P_pred * state.H' + state.R;
        K = (P_pred * state.H') / S;
        x_pred = x_pred + K * y;
        P_pred = (eye(4) - K * state.H) * P_pred;
        x_pred(3) = wrap360(x_pred(3));
        state.n_track = state.n_track + 1;
        state.n_miss  = 0;
        accepted = true;
    else
        state.n_miss = state.n_miss + 1;
    end
else
    % No detection this step. Coasting is legitimate and is why the model is
    % here, but an unbounded coast is not: after max_misses the track is
    % dropped so a stale velocity cannot steer the null indefinitely.
    state.n_miss = state.n_miss + 1;
end

if state.n_miss > cfg.max_misses
    state.x = []; state.P = [];
    state.n_track = 0; state.n_miss = 0; state.valid = false;
    return
end

state.x = x_pred;
state.P = P_pred;
state.valid = state.n_track >= cfg.min_track_steps;
pred.accepted = accepted;

% ────────────────────────── LEAD EXTRAPOLATION ────────────────────
x_lead = state.x;
for i = 1:cfg.lead_steps
    x_lead = state.F * x_lead;
end

speed = hypot(state.x(2), state.x(4) * sind(max(min(state.x(1), 180), 0)));
pred.speed_deg_s   = speed;
pred.theta_now_deg = state.x(1);
pred.phi_now_deg   = wrap360(state.x(3));
pred.theta_deg     = x_lead(1);
pred.phi_deg       = wrap360(x_lead(3));

% The speed gate. Below it the jammer is effectively static, the reactive
% covariance path is already at the oracle, and there is nothing to win by
% overriding it — so the caller is told to keep its existing branch.
pred.valid = state.valid && speed >= cfg.min_speed_deg_s && ...
    isfinite(pred.theta_deg) && isfinite(pred.phi_deg);

% A predicted elevation outside [0, 180] is off the physical grid: reflect it
% back. This is the same fold the array itself applies, so it is the right
% continuation rather than a clamp, and it keeps the lead usable when the
% jammer drifts through a pole.
if pred.theta_deg < 0 || pred.theta_deg > 180
    [pred.theta_deg, pred.phi_deg] = reflect_pole(pred.theta_deg, pred.phi_deg);
end
end


% ────────────────────────── HELPERS ───────────────────────────────

function [z, used_mirror] = fold_measurement(state, x_pred, th_m, ph_m)
% Pick the measurement branch nearer the prediction, on mirror-degenerate
% arrays only. Where the fold applies the two branches have IDENTICAL steering
% vectors, so choosing between them cannot change the null that results — it
% only keeps the track continuous, which is what the velocity estimate needs.
z = [th_m; wrap360(ph_m)];
used_mirror = false;
if ~state.mirror_fold
    return
end
d_direct = hypot(th_m - x_pred(1), wrap180(wrap360(ph_m) - x_pred(3)));
th_mir   = 180 - th_m;
d_mirror = hypot(th_mir - x_pred(1), wrap180(wrap360(ph_m) - x_pred(3)));
if d_mirror < d_direct
    z = [th_mir; wrap360(ph_m)];
    used_mirror = true;
end
end


function a = wrap360(a)
% Wrap an azimuth into [0, 360).
a = mod(a, 360);
end


function d = wrap180(d)
% Wrap an angular DIFFERENCE into [-180, 180).
d = mod(d + 180, 360) - 180;
end


function [th, ph] = reflect_pole(th, ph)
% Reflect an elevation that has run off [0, 180] back onto the sphere, taking
% the azimuth through the pole with it.
th = mod(th, 360);
if th > 180
    th = 360 - th;
    ph = wrap360(ph + 180);
end
end
