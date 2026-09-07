function state = adapt_cv_init(cv_config, dt_s, E1, E2, theta_deg, phi_deg)
% ADAPT_CV_INIT  State for the constant-velocity (theta, phi) jammer tracker.
%
%   state = ADAPT_CV_INIT(cv_config, dt_s, E1, E2, theta_deg, phi_deg)
%
%   Initializes the CV-Kalman DoA predictor that closes P12's drift gap. The
%   filter tracks [theta; theta_dot; phi; phi_dot] in DEGREES and deg/s from
%   the MUSIC angle estimates adapt_predict_update already produces, and
%   predicts the jammer's angle cv_lead_steps into the future so the null is
%   aimed where the jammer WILL be rather than where the covariance remembers
%   it having been.
%
%   WHY THIS EXISTS (measured, not assumed — [P12b] 2026-09-07)
%   On a drifting jammer the shipped reactive tracker is reproduced, to within
%   a point, by a truth-steered null delayed 10 steps — and 10 steps is exactly
%   1/(1-lambda), the covariance forgetting horizon. Steering the same null at
%   the TRUE CURRENT angle instead lifts the oracle-tracking score from 44.8%
%   to 94.8% on that case. So the drift deficit is an angle-LAG deficit, and an
%   angle predictor is the matching fix. Re-tuning lambda is not: the P12
%   sweep found no lambda where the failing cells pass.
%
%   THE MIRROR FOLD, and why this filter has to know about it
%   Some exported arrays are exactly degenerate under (theta, phi) ->
%   (180-theta, phi): data/ManyDipoles is a planar array in the z = 0 plane
%   with a common phase centre and no geometric phase re-added, and its
%   steering vectors for the two angles are IDENTICAL to 2.8e-6 relative. Two
%   consequences, and they pull opposite ways:
%     - benign for NULLING: the two steering columns being identical means a
%       null at the mirror angle IS the null at the true angle;
%     - fatal for TRACKING: the MUSIC peak hops between the two branches, and
%       a naive CV filter differences the two hops into a huge phantom
%       velocity, which is worse than no predictor at all.
%   So the filter measures the array's own mirror coherence ONCE at init (it
%   is a property of the geometry, not of the run) and, when the array is
%   degenerate, folds each measurement to whichever branch is nearer the
%   current prediction before using it. On a non-degenerate array the fold is
%   disabled and the measurement is taken as given.
%
%   Inputs:
%       cv_config : struct, the adapt.predict.cv block. ALL keys required
%                   (CLAUDE.md rule 4 — no silent defaults):
%                     lead_steps        : steps to predict ahead [steps].
%                     q_accel_deg_s2    : process noise, angular acceleration
%                                         std [deg/s^2].
%                     r_meas_deg        : DoA measurement noise std [deg].
%                     gate_deg          : innovation gate; a measurement
%                                         further than this from the
%                                         prediction is an outlier [deg].
%                     min_track_steps   : consecutive accepted updates before
%                                         the prediction may steer the null.
%                     min_speed_deg_s   : the prediction steers the null only
%                                         above this angular speed [deg/s].
%                                         Below it the caller keeps its
%                                         existing branch, so a STATIC or
%                                         on/off run is unaffected.
%                     max_misses        : consecutive gated-out measurements
%                                         before the track is dropped and
%                                         re-initialized.
%       dt_s      : closed-loop step [s] (sim.dt_s).
%       E1        : (N_el x N_theta x N_phi) primary far-field stack.
%       E2        : (N_el x N_theta x N_phi) secondary component, or [].
%       theta_deg : (1 x N_theta) elevation grid [deg].
%       phi_deg   : (1 x N_phi) azimuth grid [deg].
%
%   Outputs:
%       state : struct with fields
%           x, P        : (4x1) state [theta; theta_dot; phi; phi_dot] and its
%                         (4x4) covariance. Empty x until the first fix.
%           F, Q, H, R  : the CV model matrices (constant, built here).
%           valid       : true once min_track_steps have been accepted.
%           n_track     : consecutive accepted updates.
%           n_miss      : consecutive gated-out measurements.
%           mirror_fold : true if the array is mirror-degenerate (see above).
%           mirror_coh  : the measured coherence, kept for the report.
%           cfg         : the validated cv_config.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P12b].

% ────────────────────────── CONFIG VALIDATION ─────────────────────
REQ = {'lead_steps', 'q_accel_deg_s2', 'r_meas_deg', 'gate_deg', ...
       'min_track_steps', 'min_speed_deg_s', 'max_misses'};
for i = 1:numel(REQ)
    if ~isfield(cv_config, REQ{i}) || isempty(cv_config.(REQ{i}))
        error('adapt_cv_init:MissingKey', ...
            'Missing required adapt.predict.cv key: ''%s''.', REQ{i});
    end
end
if dt_s <= 0
    error('adapt_cv_init:BadStep', 'dt_s must be positive; got %g.', dt_s);
end

state = struct();
state.cfg  = cv_config;
state.dt_s = dt_s;

% ────────────────────────── CV MODEL MATRICES ─────────────────────
% Two decoupled 1-D constant-velocity chains (theta, phi). They are kept in
% ONE 4-state filter rather than two 2-state ones only so the gate can be
% applied to the joint innovation, which is what the angular gate means.
dt = dt_s;
state.F = [1 dt 0  0;
           0  1 0  0;
           0  0 1 dt;
           0  0 0  1];
% Discrete white-noise acceleration: the standard [dt^4/4, dt^3/2; dt^3/2,
% dt^2] block scaled by the acceleration variance. A drifting jammer is
% modelled as constant-velocity with this much unmodelled manoeuvre.
q  = cv_config.q_accel_deg_s2^2;
Qb = q * [dt^4/4, dt^3/2; dt^3/2, dt^2];
state.Q = blkdiag(Qb, Qb);
state.H = [1 0 0 0;
           0 0 1 0];
state.R = diag([cv_config.r_meas_deg^2, cv_config.r_meas_deg^2]);

% ────────────────────────── TRACK STATE ───────────────────────────
state.x       = [];        % no fix yet
state.P       = [];
state.valid   = false;
state.n_track = 0;
state.n_miss  = 0;

% ────────────────────────── MIRROR-FOLD DETECTION ─────────────────
% Coherence |e_a' e_b| / (|e_a||e_b|) between a handful of angles and their
% theta-mirrors. A planar array with no geometric phase re-added returns
% exactly 1; an array with a ground plane returns well below it. Probing a few
% angles is enough because the degeneracy is a symmetry of the whole export,
% not a local accident — the minimum over the probes is used so a single
% coincidental match cannot switch the fold on.
probe_theta = [110, 130, 150];
probe_phi   = [0, 90, 200, 300];
coh_min = Inf;
for it = 1:numel(probe_theta)
    for ip = 1:numel(probe_phi)
        c = mirror_coherence(E1, E2, theta_deg, phi_deg, ...
            probe_theta(it), probe_phi(ip));
        if isfinite(c), coh_min = min(coh_min, c); end
    end
end
if ~isfinite(coh_min), coh_min = 0; end
state.mirror_coh  = coh_min;
% 0.99 rather than 1.0: the export carries ~1e-5 of numerical asymmetry, and
% the measured split between the arrays in data/ is unambiguous (1.0000 for
% the degenerate one against 0.24-0.87 for the other two), so the threshold
% is not near anything.
state.mirror_fold = coh_min > 0.99;
end


% ────────────────────────── HELPERS ───────────────────────────────

function c = mirror_coherence(E1, E2, theta_deg, phi_deg, th_q, ph_q)
% Coherence between the steering column at (th_q, ph_q) and its theta-mirror.
[it1, ip1] = nearest_index_2d(theta_deg, phi_deg, th_q, ph_q);
[it2, ip2] = nearest_index_2d(theta_deg, phi_deg, 180 - th_q, ph_q);
a = E1(:, it1, ip1); b = E1(:, it2, ip2);
a = a(:); b = b(:);
if ~isempty(E2)
    a2 = E2(:, it1, ip1); b2 = E2(:, it2, ip2);
    a = [a; a2(:)];
    b = [b; b2(:)];
end
na = norm(a); nb = norm(b);
if na == 0 || nb == 0
    c = NaN;
    return
end
c = abs(a' * b) / (na * nb);
end
