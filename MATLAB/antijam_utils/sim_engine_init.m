function sim_state = sim_engine_init(E_primary, E_secondary, angle_deg, ...
                                     scenario, antijam_config, sim_config, mode)
% SIM_ENGINE_INIT  Build the closed-loop simulation state struct.
%
%   sim_state = SIM_ENGINE_INIT(E_primary, E_secondary, angle_deg, scenario, ...
%                               antijam_config, sim_config, mode)
%
%   Precomputes steering data and power levels for sim_engine_step. The engine
%   operates on a 1-D azimuth cut: callers extract E_primary/E_secondary from
%   the (N_el, N_theta, N_phi) stack per config antijam.cut_type (see
%   run_antijam). Polarization convention is config-selectable and enters ONLY
%   through whether E_secondary is empty (copol / single component) or not
%   ('total': incoherent sum over the two components) — mirroring
%   run_optimizer's element_patterns_secondary argument.
%
%   Power reference: sigma_n^2 = 1 (0 dB noise floor per element);
%   sigma_j^2 = 10^(jn_ratio_db/10); sigma_s^2 = 10^(sigma_s_db/10).
%   With a secondary component, signal and jammer power are split equally
%   between the two components (unpolarized-source assumption, independent
%   complex-Gaussian amplitudes per component).
%
%   The true jammer angle is mapped to the NEAREST column of the cut grid
%   (nearest_index), so angular resolution is limited by the CST export step.
%
%   Inputs:
%       E_primary      : (N_el x N_angles) complex cut of the primary component.
%       E_secondary    : (N_el x N_angles) complex cut of the secondary
%                        component, or [] for single-component operation.
%       angle_deg      : (1 x N_angles) cut angle axis [deg].
%       scenario       : struct from sim_scenario.
%       antijam_config : struct. Required: theta_s_deg, sigma_s_db.
%       sim_config     : struct. Required: seed; snapshots_per_step (Mode C).
%       mode           : 'C' (obs carries snapshots) | 'S' (scalar SINR only).
%
%   Outputs:
%       sim_state : struct with fields
%           E1, E2      : the cut matrices ([] E2 when single-component).
%           angle_deg   : (1 x N_angles) cut axis.
%           e_s         : (N_el x n_comp) steering column(s) at theta_s_deg.
%           scenario    : as passed in.
%           mode        : 'C' | 'S'.
%           k           : current step index (0 = not stepped yet; advanced by
%                         sim_engine_step, 1-based into scenario arrays).
%           sigma_s_sq, sigma_n_sq : linear powers (jammer power is per-step
%                         from scenario.jn_ratio_db).
%           n_snapshots : snapshots per step (Mode C; 0 in Mode S).
%           stream      : RandStream('mt19937ar', seed) — private stream so
%                         Monte Carlo runs are independent of the global rng.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P1].

req_field(antijam_config, 'theta_s_deg', 'antijam');
req_field(antijam_config, 'sigma_s_db',  'antijam');
req_field(sim_config, 'seed', 'sim');

if ~(ischar(mode) && any(strcmp(mode, {'C', 'S'})))
    error('sim_engine_init:BadMode', 'mode must be ''C'' or ''S''.');
end
if strcmp(mode, 'C')
    n_snapshots = req_field(sim_config, 'snapshots_per_step', 'sim');
else
    n_snapshots = 0;
end

angle_deg = angle_deg(:).';
n_angles  = numel(angle_deg);
if size(E_primary, 2) ~= n_angles
    error('sim_engine_init:BadShape', ...
        'E_primary has %d columns but angle_deg has %d entries.', ...
        size(E_primary, 2), n_angles);
end
if ~isempty(E_secondary) && ~isequal(size(E_secondary), size(E_primary))
    error('sim_engine_init:BadShape', ...
        'E_secondary size must match E_primary or be [].');
end

REQUIRED_SCENARIO = {'t_s', 'theta_j_deg', 'jammer_on', 'jn_ratio_db'};
for i = 1:numel(REQUIRED_SCENARIO)
    if ~isfield(scenario, REQUIRED_SCENARIO{i})
        error('sim_engine_init:BadScenario', ...
            'scenario struct is missing field ''%s''.', REQUIRED_SCENARIO{i});
    end
end

idx_s = nearest_index(angle_deg(:), antijam_config.theta_s_deg);
e_s   = E_primary(:, idx_s);
if ~isempty(E_secondary)
    e_s = [e_s, E_secondary(:, idx_s)];
end

sim_state = struct( ...
    'E1',          E_primary, ...
    'E2',          E_secondary, ...
    'angle_deg',   angle_deg, ...
    'e_s',         e_s, ...
    'scenario',    scenario, ...
    'mode',        mode, ...
    'k',           0, ...
    'sigma_s_sq',  10^(antijam_config.sigma_s_db / 10), ...
    'sigma_n_sq',  1.0, ...
    'n_snapshots', n_snapshots, ...
    'stream',      RandStream('mt19937ar', 'Seed', sim_config.seed));
end


function v = req_field(s, key, section)
% Fetch a required key or raise a descriptive error (no silent defaults).
if ~isfield(s, key) || isempty(s.(key))
    error('sim_engine_init:MissingKey', ...
        'Missing required %s config key: ''%s''.', section, key);
end
v = s.(key);
end
