function [jammer_state, history] = detect_jammer(covariance, history, array, ...
                                                 covariance_horizon_steps)
% ══════════════════════════════════════════════════════════════════
% DETECT_JAMMER
% Where is it, and what is it doing? (the "detect" half of Approach 1)
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   [jammer_state, history] = DETECT_JAMMER(covariance, history, array, ...
%                                           covariance_horizon_steps)
%
%   This function contains almost no logic, and that is its purpose. It asks
%   estimate_jammer_angle WHERE the jammer is at this instant, appends the
%   answer to the running history, asks classify_jammer_motion WHAT the jammer
%   has been doing, and returns both as one struct.
%
%   The whole of Approach 1's detection stage therefore reads, in the main
%   script, as a single line whose meaning is obvious -- and anyone who wants
%   the detail opens one of two small files, each answering one question with
%   one kind of evidence. The previous implementation fused both questions plus
%   the weight solve into one 485-line function; this seam is what prevents that.
%
%   IT NEVER SEES THE TRUTH. Its inputs are a covariance estimated from
%   snapshots, its own past outputs, and the receiver's knowledge of its own
%   antenna. The scenario struct is not among them, and must never be.
%
%   Inputs:
%       covariance               : (n_elements x n_elements) from
%                                  sample_covariance. Units: power.
%       history                  : struct of past estimates, or [] on the first
%                                  step. Returned updated -- pass it back in.
%       array                    : struct from make_array.
%       covariance_horizon_steps : 1/(1 - forgetting_lambda). Units: steps.
%
%   Outputs:
%       jammer_state : struct with fields
%           is_feasible   : logical, the array can do direction finding at all.
%           is_present    : logical, a jammer was detected this step.
%           behaviour     : 'steady' | 'onoff' | 'drifting' | 'absent'.
%           theta_deg     : where to place the null. Units: degrees.
%           phi_deg       : ditto. Units: degrees.
%           measured_theta_deg : this step's raw beamscan estimate, before any
%                                averaging or extrapolation. Units: degrees.
%           music_theta_deg    : the same step's MUSIC estimate. Reported for
%                                comparison only -- it never places a null.
%                                Units: degrees.
%           drift_rate_deg_step : fitted angular rate. Units: deg/step.
%           lead_deg      : extrapolation applied, 0 unless drifting. Units: deg.
%           duty_cycle    : fraction of the recent window the jammer was seen.
%           peak_ratio_db : MUSIC peak sharpness this step. Units: dB.
%           reason        : one line explaining the verdict, for reporting.
%       history      : updated history struct; feed it back on the next step.

% ────────────────────────── WHERE IS IT NOW? ──────────────────────

jammer_estimate = estimate_jammer_angle(covariance, array);

% ────────────────────────── REMEMBER IT ───────────────────────────

if isempty(history)
    history = struct('theta_deg', [], 'phi_deg', [], 'is_present', []);
end

history.theta_deg  = [history.theta_deg;  jammer_estimate.theta_deg];
history.phi_deg    = [history.phi_deg;    jammer_estimate.phi_deg];
history.is_present = [history.is_present; jammer_estimate.is_present];

% ────────────────────────── WHAT HAS IT BEEN DOING? ───────────────

motion = classify_jammer_motion(history.theta_deg, history.phi_deg, ...
                                history.is_present, array.profile, ...
                                covariance_horizon_steps);

% ────────────────────────── ASSEMBLE ──────────────────────────────

% is_present here is the CLASSIFIER's verdict, not this step's raw detection: a
% jammer that is switched off right now is still a jammer to be nulled, because
% it will come back and the null costs almost nothing while it is away.
jammer_state = struct( ...
    'is_feasible',         jammer_estimate.is_feasible, ...
    'is_present',          ~strcmp(motion.behaviour, 'absent'), ...
    'behaviour',           motion.behaviour, ...
    'theta_deg',           motion.aim_theta_deg, ...
    'phi_deg',             motion.aim_phi_deg, ...
    'measured_theta_deg',  jammer_estimate.theta_deg, ...
    'music_theta_deg',     jammer_estimate.music_theta_deg, ...
    'drift_rate_deg_step', motion.drift_rate_deg_step, ...
    'lead_deg',            motion.lead_deg, ...
    'duty_cycle',          motion.duty_cycle, ...
    'peak_ratio_db',       jammer_estimate.peak_ratio_db, ...
    'reason',              motion.reason);

% An array that cannot do direction finding at all reports that, rather than
% the classifier's opinion about a history of NaNs.
if ~jammer_estimate.is_feasible
    jammer_state.is_present = false;
    jammer_state.behaviour  = 'absent';
    jammer_state.reason     = jammer_estimate.reason;
end
end
