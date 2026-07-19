function codebook = agent_codebook_build(element_patterns_stacked, ...
    element_patterns_secondary, theta_deg, phi_deg, antijam_config, ...
    agent_config, optimizer_config, cache_path)
% AGENT_CODEBOOK_BUILD  Synthesize (or load cached) codebook of nulling arms.
%
%   codebook = AGENT_CODEBOOK_BUILD(element_patterns_stacked, ...
%       element_patterns_secondary, theta_deg, phi_deg, antijam_config, ...
%       agent_config, optimizer_config, cache_path)
%
%   One arm per candidate null sector: calls the Milestone-1 run_optimizer with
%   two directives — peak at theta_s_deg (width agent_config.peak_width_deg,
%   weight 1) + null sector centered on the arm's angle (width
%   agent_config.null_width_deg, weight agent_config.null_weight) — on the
%   FULL 2-D grid (the same directives convention as run_optimization), then
%   POLISHES each arm by null-space projection: w <- (I - U*U') * w with U an
%   orthonormal basis of the null-window steering columns (both polarization
%   components when a secondary stack is given). The projection places exact
%   zeros at the sampled window angles; the optimizer solution provides the
%   well-shaped beam the projection perturbs only slightly. (P4 finding: the
%   composite cost alone plateaus near -25 dB null depth — a true Pareto
%   point, not a solver artifact.)
%   Arm 1 is the UNCONSTRAINED peak-only solution (null_center_deg = NaN): it
%   serves jammer-off periods and is the peak-gain reference for the P4 gate.
%   Null centers span the cut axis on an agent_config.null_grid_deg grid,
%   excluding the main-beam guard sector.
%
%   Cut-axis angles map to full-grid directive coordinates as:
%       'phi_cut'   : theta = cut_theta_deg, phi = mod(angle, 360)
%       'theta_cut' : theta = |angle|, phi = cut_phi_deg (+180 for angle < 0)
%
%   Offline and cached: if cache_path exists and its stored build parameters
%   match the current config (isequaln on the params echo), the codebook is
%   loaded instead of rebuilt. Reproducible: run_optimizer restarts are
%   internally seeded.
%
%   Inputs:
%       element_patterns_stacked   : (N_el, N_theta, N_phi) complex stack.
%       element_patterns_secondary : secondary-component stack or [] (same
%                                    polarization convention as run_optimizer).
%       theta_deg, phi_deg         : full angle grids [deg].
%       antijam_config             : struct. Required: theta_s_deg, guard_deg,
%                                    cut_type, and cut_theta_deg (phi_cut) or
%                                    cut_phi_deg (theta_cut).
%       agent_config               : struct. Required: null_grid_deg,
%                                    null_width_deg, peak_width_deg, null_weight.
%       optimizer_config           : struct, per run_optimizer requirements.
%       cache_path                 : char, .mat file for the cached codebook;
%                                    [] disables caching.
%
%   Outputs:
%       codebook : struct with fields
%           W               : (N_el x n_arms) complex weight columns.
%           null_center_deg : (1 x n_arms) arm null centers on the cut axis
%                             (NaN for the unconstrained arm 1).
%           meta            : struct with params (build-parameter echo used
%                             for cache validation) and built (timestamp).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P4].

if nargin < 8
    cache_path = [];
end

theta_s = req_field(antijam_config, 'theta_s_deg', 'antijam');
guard   = req_field(antijam_config, 'guard_deg',   'antijam');
cut_type = req_field(antijam_config, 'cut_type',   'antijam');
null_grid  = req_field(agent_config, 'null_grid_deg',  'agent');
null_width = req_field(agent_config, 'null_width_deg', 'agent');
peak_width = req_field(agent_config, 'peak_width_deg', 'agent');
null_weight = req_field(agent_config, 'null_weight',   'agent');

% ── Candidate null centers on the cut axis, guard sector excluded ─
switch cut_type
    case 'phi_cut'
        cut_fixed = req_field(antijam_config, 'cut_theta_deg', 'antijam');
        cand = 0:null_grid:(360 - null_grid / 2);
        dist = abs(mod(cand - theta_s + 180, 360) - 180);   % circular distance
    case 'theta_cut'
        cut_fixed = req_field(antijam_config, 'cut_phi_deg', 'antijam');
        cand = -90:null_grid:90;
        dist = abs(cand - theta_s);
    otherwise
        error('agent_codebook_build:BadCutType', ...
            'Unknown cut_type ''%s'' (expected ''phi_cut'' or ''theta_cut'').', cut_type);
end
centers = cand(dist >= guard);
n_arms  = numel(centers) + 1;

% A null window overlapping the peak window makes the two directives fight;
% the optimizer then sacrifices one of them (observed in P4 tuning).
min_guard = (peak_width + null_width) / 2;
if guard < min_guard
    warning('agent_codebook_build:GuardTooSmall', ...
        ['guard_deg = %g < (peak_width_deg + null_width_deg)/2 = %g: ' ...
         'edge arms overlap the peak window and will be shallow.'], ...
        guard, min_guard);
end

% ── Cache check ───────────────────────────────────────────────────
params = struct( ...
    'theta_s_deg', theta_s, 'guard_deg', guard, 'cut_type', cut_type, ...
    'cut_fixed_deg', cut_fixed, 'null_grid_deg', null_grid, ...
    'null_width_deg', null_width, 'peak_width_deg', peak_width, ...
    'null_weight', null_weight, 'optimizer', optimizer_config, ...
    'n_elements', size(element_patterns_stacked, 1), ...
    'has_secondary', ~isempty(element_patterns_secondary), ...
    'grid_size', [numel(theta_deg), numel(phi_deg)]);

if ~isempty(cache_path) && isfile(cache_path)
    cached = load(cache_path);
    if isfield(cached, 'codebook') && isequaln(cached.codebook.meta.params, params)
        codebook = cached.codebook;
        fprintf('Codebook: loaded %d arms from cache %s\n', ...
            size(codebook.W, 2), cache_path);
        return;
    end
    fprintf('Codebook: cache %s stale (parameters changed), rebuilding...\n', cache_path);
end

% ── Synthesize arms via the Milestone-1 optimizer ─────────────────
% Peak uses aggregation 'min' (maximize the WORST point in the window): with
% 'mean', the solid-angle weighting lets the optimizer park the beam at the
% window edge, off theta_s entirely (P4 tuning finding). Consequence:
% peak_width_deg must be ~ the achievable beamwidth, or 'min' forces an
% unrealizably fat beam and sacrifices gain.
[th_s, ph_s] = cut_to_theta_phi(theta_s, cut_type, cut_fixed);
peak_dir = struct('type', 'peak', 'theta', th_s, 'phi', ph_s, ...
                  'width', peak_width, 'weight', 1.0, 'aggregation', 'min');

W = complex(zeros(size(element_patterns_stacked, 1), n_arms));
null_center_deg = NaN(1, n_arms);

fprintf('Codebook: building %d arms (1 unconstrained + %d null sectors)...\n', ...
    n_arms, numel(centers));
result = run_optimizer(element_patterns_stacked, theta_deg, phi_deg, ...
    {peak_dir}, optimizer_config, element_patterns_secondary);
W(:, 1) = result.weights_complex;

for i = 1:numel(centers)
    [th_n, ph_n] = cut_to_theta_phi(centers(i), cut_type, cut_fixed);
    null_dir = struct('type', 'null', 'theta', th_n, 'phi', ph_n, ...
                      'width', null_width, 'weight', null_weight);
    result = run_optimizer(element_patterns_stacked, theta_deg, phi_deg, ...
        {peak_dir, null_dir}, optimizer_config, element_patterns_secondary);
    w_arm = result.weights_complex;

    % Null-space projection polish: exact zeros on the sampled window angles.
    % The mask is restricted to the 1-D CUT the engine operates on (the
    % milestone is 2-D azimuth): projecting the full 2-D theta x phi window
    % nulls angles the simulation never uses and eats DOF quadratically —
    % on the real 5-deg grid a 20-deg window spanned rank 15 of 20 DOF and
    % destroyed the beam (P6 finding).
    mask = angular_window_mask(theta_deg(:), phi_deg(:), th_n, ph_n, ...
                               null_width, null_width);
    mask = restrict_mask_to_cut(mask, theta_deg, phi_deg, cut_type, cut_fixed);
    E_win = window_steering_columns(element_patterns_stacked, mask);
    if ~isempty(element_patterns_secondary)
        E_win = [E_win, window_steering_columns(element_patterns_secondary, mask)]; %#ok<AGROW>
    end
    [U, S, ~] = svd(E_win, 'econ');
    sv = diag(S);
    r  = sum(sv > 1e-6 * sv(1));
    if r > size(W, 1) / 2
        warning('agent_codebook_build:WindowRankHigh', ...
            'Null window at %g deg spans rank %d of %d DOF; beam quality will suffer.', ...
            centers(i), r, size(W, 1));
    end
    w_arm = w_arm - U(:, 1:r) * (U(:, 1:r)' * w_arm);

    W(:, i + 1) = w_arm;
    null_center_deg(i + 1) = centers(i);
end

codebook = struct('W', W, 'null_center_deg', null_center_deg, ...
    'meta', struct('params', params, ...
                   'built', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'))));

if ~isempty(cache_path)
    cache_dir = fileparts(cache_path);
    if ~isempty(cache_dir) && ~isfolder(cache_dir)
        mkdir(cache_dir);
    end
    save(cache_path, 'codebook');
    fprintf('Codebook: cached to %s\n', cache_path);
end
end


% ────────────────────────── HELPERS ───────────────────────────────

function v = req_field(s, key, section)
% Fetch a required key or raise a descriptive error (no silent defaults).
if ~isfield(s, key) || isempty(s.(key))
    error('agent_codebook_build:MissingKey', ...
        'Missing required %s config key: ''%s''.', section, key);
end
v = s.(key);
end


function mask = restrict_mask_to_cut(mask, theta_deg, phi_deg, cut_type, cut_fixed)
% Keep only the window points that lie ON the engine's 1-D cut.
if strcmp(cut_type, 'phi_cut')
    it = nearest_index(theta_deg(:), cut_fixed);
    keep = false(size(mask));
    keep(it, :) = true;
else
    ip0   = nearest_index(phi_deg(:), mod(cut_fixed, 360));
    ip180 = nearest_index(phi_deg(:), mod(cut_fixed + 180, 360));
    keep = false(size(mask));
    keep(:, [ip0, ip180]) = true;
end
mask = mask & keep;
end


function E_win = window_steering_columns(stack, mask)
% Collect the (N_el x M) steering columns of the M grid points where the
% (N_theta x N_phi) window mask is true.
n_el  = size(stack, 1);
flat  = reshape(stack, n_el, []);        % (N_el x N_theta*N_phi), grid flattened
E_win = flat(:, mask(:));
end


function [th, ph] = cut_to_theta_phi(angle_deg, cut_type, cut_fixed)
% Map a cut-axis angle to (theta, phi) directive coordinates on the full grid.
if strcmp(cut_type, 'phi_cut')
    th = cut_fixed;
    ph = mod(angle_deg, 360);
else
    th = abs(angle_deg);
    ph = mod(cut_fixed + 180 * (angle_deg < 0), 360);
end
end
