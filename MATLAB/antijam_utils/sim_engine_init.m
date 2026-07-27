function sim_state = sim_engine_init(stack1, stack2, theta_deg, phi_deg, ...
                                     scenario, antijam_config, sim_config, mode)
% SIM_ENGINE_INIT  Build the closed-loop simulation state struct.
%
%   sim_state = SIM_ENGINE_INIT(stack1, stack2, theta_deg, phi_deg, scenario, ...
%                               antijam_config, sim_config, mode)
%
%   Precomputes steering data and power levels for sim_engine_step. [P7]: the
%   engine operates natively on the full 2-D (theta, phi) far-field grid — no
%   1-D cut is extracted. Polarization convention is config-selectable and
%   enters ONLY through whether stack2 is empty (copol / single component) or
%   not ('total': incoherent sum over the two components) — mirroring
%   run_optimizer's element_patterns_secondary argument.
%
%   Power reference: sigma_n^2 = 1 (0 dB noise floor per element);
%   sigma_j^2 = 10^(jn_ratio_db/10); sigma_s^2 = 10^(sigma_s_db/10).
%   With a secondary component, signal and jammer power are split equally
%   between the two components (unpolarized-source assumption, independent
%   complex-Gaussian amplitudes per component).
%
%   The true jammer (theta, phi) is mapped to the NEAREST grid point
%   (nearest_index_2d), so angular resolution is limited by the CST export
%   step. The desired direction (theta_s_deg, phi_s_deg) likewise maps to a
%   single nearest grid point (matching the pre-[P7] theta_s_deg-only behavior,
%   which also had no explicit peak width).
%
%   Inputs:
%       stack1, stack2 : (N_el x N_theta x N_phi) complex far-field stacks
%                        (primary / secondary polarization component); stack2
%                        may be [] for single-component operation.
%       theta_deg      : (1 x N_theta) elevation grid [deg].
%       phi_deg        : (1 x N_phi) azimuth grid [deg].
%       scenario       : struct from sim_scenario.
%       antijam_config : struct. Required: theta_s_deg, phi_s_deg, sigma_s_db.
%       sim_config     : struct. Required: seed; snapshots_per_step (Mode C).
%       mode           : 'C' (obs carries snapshots) | 'S' (scalar SINR only).
%
%   Outputs:
%       sim_state : struct with fields
%           E1, E2      : (N_el x N_theta*N_phi) flattened stacks ([] E2 when
%                         single-component), column-major over (theta, phi)
%                         i.e. matching reshape(stack, N_el, N_theta*N_phi).
%           theta_deg, phi_deg : grid axes.
%           e_s         : (N_el x n_comp) steering column(s) at
%                         (theta_s_deg, phi_s_deg).
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
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P1, P7].

req_field(antijam_config, 'theta_s_deg', 'antijam');
req_field(antijam_config, 'phi_s_deg',   'antijam');
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

theta_deg = theta_deg(:).';
phi_deg   = phi_deg(:).';
n_theta   = numel(theta_deg);
n_phi     = numel(phi_deg);
n_el      = size(stack1, 1);
% [size(stack1,3) rather than size(stack1) avoids MATLAB silently dropping a
% trailing singleton dimension (n_phi = 1 degenerate grids, e.g. toy tests).]
if size(stack1, 1) ~= n_el || size(stack1, 2) ~= n_theta || size(stack1, 3) ~= n_phi
    error('sim_engine_init:BadShape', ...
        'stack1 must be (N_el x %d x %d); got %s.', ...
        n_theta, n_phi, mat2str(size(stack1)));
end
if ~isempty(stack2) && ...
        (size(stack2, 1) ~= n_el || size(stack2, 2) ~= n_theta || size(stack2, 3) ~= n_phi)
    error('sim_engine_init:BadShape', ...
        'stack2 size must match stack1 or be [].');
end

REQUIRED_SCENARIO = {'t_s', 'theta_j_deg', 'phi_j_deg', 'jammer_on', 'jn_ratio_db'};
for i = 1:numel(REQUIRED_SCENARIO)
    if ~isfield(scenario, REQUIRED_SCENARIO{i})
        error('sim_engine_init:BadScenario', ...
            'scenario struct is missing field ''%s''.', REQUIRED_SCENARIO{i});
    end
end

E1 = reshape(stack1, n_el, n_theta * n_phi);
E2 = [];
if ~isempty(stack2)
    E2 = reshape(stack2, n_el, n_theta * n_phi);
end

[it_s, ip_s] = nearest_index_2d(theta_deg, phi_deg, ...
    antijam_config.theta_s_deg, antijam_config.phi_s_deg);
idx_s = (ip_s - 1) * n_theta + it_s;   % linear index into the (theta,phi) flatten
e_s   = E1(:, idx_s);
if ~isempty(E2)
    e_s = [e_s, E2(:, idx_s)];
end

sim_state = struct( ...
    'E1',          E1, ...
    'E2',          E2, ...
    'theta_deg',   theta_deg, ...
    'phi_deg',     phi_deg, ...
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
