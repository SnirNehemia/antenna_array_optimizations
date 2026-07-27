function scenario = sim_scenario(scenario_config, antijam_config, sim_config)
% SIM_SCENARIO  Generate a deterministic jammer trajectory + power timeline.
%
%   scenario = SIM_SCENARIO(scenario_config, antijam_config, sim_config)
%
%   Expands one scenario definition (config antijam.scenarios) into a fully
%   sampled timeline on the closed-loop time grid. Deterministic under
%   sim_config.seed: the static/initial jammer position is drawn per seed,
%   uniform over the full sphere EXCLUDING the main-beam guard cap (spherical
%   angular separation < antijam_config.guard_deg from (theta_s_deg, phi_s_deg),
%   via angular_separation_deg).
%
%   Motion model (native 2D, [P7]): theta and phi each move on an unreflected
%   linear path (u_theta(t) = theta0 + theta_drift_deg_per_s * t, similarly for
%   phi; optional jump once t >= jump_time_s). theta is then folded
%   (triangle-wave reflection) into the physical elevation range [0, 180]; phi
%   wraps circularly mod 360 (no folding needed). If the folded point still
%   falls inside the guard cap, it is pushed radially outward — in a local
%   tangent-plane approximation scaled by sin(theta_s) to account for azimuth
%   convergence near the pole — to sit exactly on the guard boundary. This is a
%   heuristic (not a true 2-D billiard reflection) but keeps the trajectory
%   continuous and guard-respecting, which is all the closed-loop KPIs need.
%
%   Inputs:
%       scenario_config : struct, one element of config antijam.scenarios.
%                         Required: id, motion ('static'|'drift'), power
%                         ('constant'|'step'|'onoff'|'window'). Motion 'drift'
%                         requires theta_drift_deg_per_s + phi_drift_deg_per_s;
%                         optional theta_jump_deg + phi_jump_deg + jump_time_s.
%                         Power 'step' requires power_step_db + step_time_s;
%                         'onoff' requires duty_cycle + toggle_period_s;
%                         'window' requires on_time_s + off_time_s (jammer
%                         silent outside [on, off) — the null-lifecycle case).
%                         Optional duration_s overrides sim_config.duration_s
%                         for this scenario. Optional theta_j_deg + phi_j_deg
%                         fix the static/initial jammer position instead of the
%                         per-seed random draw (error if inside the guard cap).
%                         Optional jn_ratio_db overrides antijam_config.jn_ratio_db.
%       antijam_config  : struct. Required: theta_s_deg, phi_s_deg, guard_deg,
%                         jn_ratio_db.
%       sim_config      : struct. Required: dt_s, duration_s, seed.
%
%   Outputs:
%       scenario : struct with fields
%           id          : char, scenario id (e.g. 'S2').
%           t_s         : (1 x T) time grid [s], T = floor(duration_s/dt_s) + 1.
%           theta_j_deg : (1 x T) true jammer elevation [deg].
%           phi_j_deg   : (1 x T) true jammer azimuth [deg].
%           jammer_on   : (1 x T) logical, false during 'onoff' off-phases.
%           jn_ratio_db : (1 x T) jammer-to-noise ratio [dB] (steps applied).
%           events      : cell array of structs {t_s, type} for every jammer
%                         event ('jump' | 'power_step' | 'turn_on') — consumed
%                         by kpi_evaluate for recovery-time measurement.
%
%   NOTE: theta_j_deg/phi_j_deg are ground truth for the simulator, oracle, and
%   KPI code ONLY. They must never be passed to adapt_* / agent_* algorithms.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P1, P7].

req_field(scenario_config, 'id',    'scenario');
req_field(scenario_config, 'motion', 'scenario');
req_field(scenario_config, 'power', 'scenario');
req_field(antijam_config, 'theta_s_deg', 'antijam');
req_field(antijam_config, 'phi_s_deg',   'antijam');
req_field(antijam_config, 'guard_deg',   'antijam');
req_field(antijam_config, 'jn_ratio_db', 'antijam');
req_field(sim_config, 'dt_s',       'sim');
req_field(sim_config, 'duration_s', 'sim');
req_field(sim_config, 'seed',       'sim');

dt = sim_config.dt_s;
duration_s = sim_config.duration_s;
if isfield(scenario_config, 'duration_s') && ~isempty(scenario_config.duration_s)
    duration_s = scenario_config.duration_s;   % per-scenario override (e.g. S6)
end
n_steps = floor(duration_s / dt) + 1;
t_s     = (0:n_steps - 1) * dt;

theta_s = antijam_config.theta_s_deg;
phi_s   = antijam_config.phi_s_deg;
guard   = antijam_config.guard_deg;
stream  = RandStream('mt19937ar', 'Seed', sim_config.seed);

% ── Initial position: fixed (config) or per-seed random draw outside guard ─
have_fixed = isfield(scenario_config, 'theta_j_deg') && ~isempty(scenario_config.theta_j_deg);
if have_fixed
    req_field(scenario_config, 'phi_j_deg', 'scenario');
    theta0 = scenario_config.theta_j_deg;
    phi0   = scenario_config.phi_j_deg;
    if angular_separation_deg(theta0, phi0, theta_s, phi_s) < guard
        error('sim_scenario:AngleInGuard', ...
            'Scenario %s: (theta_j_deg, phi_j_deg) = (%g, %g) is inside the guard cap (%g deg of (theta_s, phi_s)).', ...
            scenario_config.id, theta0, phi0, guard);
    end
else
    max_tries = 10000;
    theta0 = NaN; phi0 = NaN;
    for i = 1:max_tries
        cand_theta = 180 * rand(stream);
        cand_phi   = 360 * rand(stream);
        if angular_separation_deg(cand_theta, cand_phi, theta_s, phi_s) >= guard
            theta0 = cand_theta;
            phi0   = cand_phi;
            break
        end
    end
    if isnan(theta0)
        error('sim_scenario:NoJammerSpan', ...
            'Guard cap (guard_deg = %g) leaves no room for a random jammer draw after %d tries.', ...
            guard, max_tries);
    end
end

% ── Motion: unreflected linear path per axis, then folded/guard-clamped ────
events = {};
switch scenario_config.motion
    case 'static'
        theta_rate = 0.0;
        phi_rate   = 0.0;
    case 'drift'
        theta_rate = req_field(scenario_config, 'theta_drift_deg_per_s', 'scenario');
        phi_rate   = req_field(scenario_config, 'phi_drift_deg_per_s',   'scenario');
    otherwise
        error('sim_scenario:BadMotion', ...
            'Unknown motion ''%s'' (expected ''static'' or ''drift'').', ...
            scenario_config.motion);
end
raw_theta = theta0 + theta_rate * t_s;
raw_phi   = phi0   + phi_rate   * t_s;
if isfield(scenario_config, 'theta_jump_deg') && ~isempty(scenario_config.theta_jump_deg)
    jump_t = req_field(scenario_config, 'jump_time_s', 'scenario');
    mask   = t_s >= jump_t;
    raw_theta(mask) = raw_theta(mask) + scenario_config.theta_jump_deg;
    if isfield(scenario_config, 'phi_jump_deg') && ~isempty(scenario_config.phi_jump_deg)
        raw_phi(mask) = raw_phi(mask) + scenario_config.phi_jump_deg;
    end
    events{end + 1} = struct('t_s', jump_t, 'type', 'jump');
end

theta_f = reflect_into(raw_theta, 0, 180);
phi_f   = mod(raw_phi, 360);
[theta_j, phi_j] = guard_clamp(theta_f, phi_f, theta_s, phi_s, guard);

% ── Power profile ─────────────────────────────────────────────────
base_jn = antijam_config.jn_ratio_db;
if isfield(scenario_config, 'jn_ratio_db') && ~isempty(scenario_config.jn_ratio_db)
    base_jn = scenario_config.jn_ratio_db;      % per-scenario override
end
jn = base_jn * ones(1, n_steps);
on = true(1, n_steps);
switch scenario_config.power
    case 'constant'
        % nothing to do
    case 'step'
        step_db = req_field(scenario_config, 'power_step_db', 'scenario');
        step_t  = req_field(scenario_config, 'step_time_s',   'scenario');
        mask     = t_s >= step_t;
        jn(mask) = jn(mask) + step_db;
        events{end + 1} = struct('t_s', step_t, 'type', 'power_step');
    case 'onoff'
        duty   = req_field(scenario_config, 'duty_cycle',      'scenario');
        period = req_field(scenario_config, 'toggle_period_s', 'scenario');
        on = mod(t_s, period) < duty * period;
        % Every off -> on transition is a recovery event (run start excluded).
        turn_on = find(on & [false, ~on(1:end - 1)]);
        for i = 1:numel(turn_on)
            events{end + 1} = struct('t_s', t_s(turn_on(i)), 'type', 'turn_on'); %#ok<AGROW>
        end
    case 'window'
        on_t  = req_field(scenario_config, 'on_time_s',  'scenario');
        off_t = req_field(scenario_config, 'off_time_s', 'scenario');
        on = (t_s >= on_t) & (t_s < off_t);
        events{end + 1} = struct('t_s', on_t,  'type', 'turn_on');
        events{end + 1} = struct('t_s', off_t, 'type', 'turn_off');
    otherwise
        error('sim_scenario:BadPower', ...
            'Unknown power ''%s'' (expected ''constant'', ''step'', ''onoff'' or ''window'').', ...
            scenario_config.power);
end

scenario = struct( ...
    'id',          scenario_config.id, ...
    't_s',         t_s, ...
    'theta_j_deg', theta_j, ...
    'phi_j_deg',   phi_j, ...
    'jammer_on',   on, ...
    'jn_ratio_db', jn, ...
    'events',      {events});
end


% ────────────────────────── HELPERS ───────────────────────────────

function v = req_field(s, key, section)
% Fetch a required key or raise a descriptive error (no silent defaults).
if ~isfield(s, key) || isempty(s.(key))
    error('sim_scenario:MissingKey', ...
        'Missing required %s config key: ''%s''.', section, key);
end
v = s.(key);
end


function p = reflect_into(p, lo, hi)
% Fold positions into [lo, hi] by triangle-wave (billiard) reflection.
span = hi - lo;
q = mod(p - lo, 2 * span);
fold = q > span;
q(fold) = 2 * span - q(fold);
p = lo + q;
end


function [theta_out, phi_out] = guard_clamp(theta_in, phi_in, theta_s, phi_s, guard_deg)
% Push any (theta, phi) sample inside the guard cap radially out to its
% boundary, using a local tangent-plane approximation (phi scaled by
% sin(theta_s) so azimuth degrees and elevation degrees are comparable near
% the pole). Samples already outside the cap are left untouched.
theta_out = theta_in;
phi_out   = phi_in;
sep = angular_separation_deg(theta_in, phi_in, theta_s, phi_s);
viol = sep < guard_deg;
if ~any(viol)
    return
end
scale = max(sind(theta_s), 0.05);   % avoid a degenerate metric at the pole
dth = theta_in(viol) - theta_s;
dph = mod(phi_in(viol) - phi_s + 180, 360) - 180;   % shortest signed azimuth diff
dph_local = dph .* scale;
r = hypot(dth, dph_local);
degenerate = r < 1e-9;
r(degenerate) = 1;
dth(degenerate) = 1;              % arbitrary push direction if exactly on target
dph_local(degenerate) = 0;
factor = guard_deg ./ r;
dth_new       = dth .* factor;
dph_local_new = dph_local .* factor;
theta_out(viol) = theta_s + dth_new;
phi_out(viol)   = mod(phi_s + dph_local_new ./ scale, 360);
end
