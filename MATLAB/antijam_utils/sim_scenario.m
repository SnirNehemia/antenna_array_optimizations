function scenario = sim_scenario(scenario_config, antijam_config, sim_config)
% SIM_SCENARIO  Generate a deterministic jammer trajectory + power timeline.
%
%   scenario = SIM_SCENARIO(scenario_config, antijam_config, sim_config)
%
%   Expands one scenario definition (S1..S5 in config antijam.scenarios) into a
%   fully sampled timeline on the closed-loop time grid. Deterministic under
%   sim_config.seed: the static/initial jammer angle is drawn per seed, uniform
%   over the cut span EXCLUDING the main-beam guard sector
%   [theta_s_deg - guard_deg, theta_s_deg + guard_deg].
%
%   Motion model: the jammer moves on an unreflected linear path
%   u(t) = u0 + drift_deg_per_s * t (+ jump_deg once t >= jump_time_s), in the
%   offset coordinate u = theta - theta_s_deg, then the path is FOLDED
%   (triangle-wave reflection) into the allowed interval so it never enters the
%   guard sector or leaves the cut span. For cut_type 'phi_cut' the allowed
%   interval is [guard, 360 - guard] (circular span); for 'theta_cut' it is the
%   sub-interval of [-90, 90] on whichever side of the main beam the seeded
%   initial draw lands (drift stays on that side by reflection).
%
%   Inputs:
%       scenario_config : struct, one element of config antijam.scenarios.
%                         Required: id, motion ('static'|'drift'), power
%                         ('constant'|'step'|'onoff'|'window'). Motion 'drift'
%                         requires drift_deg_per_s; optional jump_deg +
%                         jump_time_s. Power 'step' requires power_step_db +
%                         step_time_s; 'onoff' requires duty_cycle +
%                         toggle_period_s; 'window' requires on_time_s +
%                         off_time_s (jammer silent outside [on, off) — the
%                         null-lifecycle case). Optional duration_s overrides
%                         sim_config.duration_s for this scenario.
%       antijam_config  : struct. Required: theta_s_deg, guard_deg,
%                         jn_ratio_db, cut_type ('phi_cut'|'theta_cut').
%       sim_config      : struct. Required: dt_s, duration_s, seed.
%
%   Outputs:
%       scenario : struct with fields
%           id          : char, scenario id (e.g. 'S2').
%           t_s         : (1 x T) time grid [s], T = floor(duration_s/dt_s) + 1.
%           theta_j_deg : (1 x T) true jammer angle along the cut axis [deg].
%           jammer_on   : (1 x T) logical, false during 'onoff' off-phases.
%           jn_ratio_db : (1 x T) jammer-to-noise ratio [dB] (steps applied).
%           events      : cell array of structs {t_s, type} for every jammer
%                         event ('jump' | 'power_step' | 'turn_on') — consumed
%                         by kpi_evaluate for recovery-time measurement.
%
%   NOTE: theta_j_deg is ground truth for the simulator, oracle, and KPI code
%   ONLY. It must never be passed to adapt_* / agent_* algorithms.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P1].

req_field(scenario_config, 'id',    'scenario');
req_field(scenario_config, 'motion', 'scenario');
req_field(scenario_config, 'power', 'scenario');
req_field(antijam_config, 'theta_s_deg', 'antijam');
req_field(antijam_config, 'guard_deg',   'antijam');
req_field(antijam_config, 'jn_ratio_db', 'antijam');
req_field(antijam_config, 'cut_type',    'antijam');
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
guard   = antijam_config.guard_deg;
stream  = RandStream('mt19937ar', 'Seed', sim_config.seed);

% ── Allowed interval [lo, hi] in offset coords u = theta - theta_s ─
switch antijam_config.cut_type
    case 'phi_cut'
        lo = guard;
        hi = 360 - guard;
        u0 = lo + (hi - lo) * rand(stream);
    case 'theta_cut'
        lo_lin  = -90 - theta_s;                 % offset span of [-90, 90]
        hi_lin  =  90 - theta_s;
        len_neg = max(0, -guard - lo_lin);       % [lo_lin, -guard]
        len_pos = max(0, hi_lin - guard);        % [guard, hi_lin]
        if len_neg + len_pos <= 0
            error('sim_scenario:NoJammerSpan', ...
                'Guard sector (guard_deg = %g) covers the whole cut span.', guard);
        end
        pick = (len_neg + len_pos) * rand(stream);
        if pick < len_neg
            lo = lo_lin;  hi = -guard;  u0 = lo + pick;
        else
            lo = guard;   hi = hi_lin;  u0 = lo + (pick - len_neg);
        end
    otherwise
        error('sim_scenario:BadCutType', ...
            'Unknown cut_type ''%s'' (expected ''phi_cut'' or ''theta_cut'').', ...
            antijam_config.cut_type);
end

% ── Motion: unreflected linear path, folded into [lo, hi] ─────────
events = {};
switch scenario_config.motion
    case 'static'
        rate = 0.0;
    case 'drift'
        rate = req_field(scenario_config, 'drift_deg_per_s', 'scenario');
    otherwise
        error('sim_scenario:BadMotion', ...
            'Unknown motion ''%s'' (expected ''static'' or ''drift'').', ...
            scenario_config.motion);
end
raw = u0 + rate * t_s;
if isfield(scenario_config, 'jump_deg') && ~isempty(scenario_config.jump_deg)
    jump_t = req_field(scenario_config, 'jump_time_s', 'scenario');
    mask   = t_s >= jump_t;
    raw(mask) = raw(mask) + scenario_config.jump_deg;
    events{end + 1} = struct('t_s', jump_t, 'type', 'jump');
end
u = reflect_into(raw, lo, hi);

if strcmp(antijam_config.cut_type, 'phi_cut')
    theta_j = mod(theta_s + u, 360);
else
    theta_j = theta_s + u;
end

% ── Power profile ─────────────────────────────────────────────────
jn = antijam_config.jn_ratio_db * ones(1, n_steps);
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
