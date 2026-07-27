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
%   two directives — peak at (theta_s_deg, phi_s_deg) (width
%   agent_config.peak_width_deg, weight 1) + null sector centered on the arm's
%   (theta, phi) (width agent_config.null_width_deg, weight
%   agent_config.null_weight) — on the FULL 2-D grid (the same directives
%   convention as run_optimization), then POLISHES each arm by null-space
%   projection: w <- w - U*(U'*w) with U the top-r left singular vectors of
%   the null-window steering columns (both polarization components when a
%   secondary stack is given), r capped at agent_config.null_rank_cap. The
%   projection places (approximately, for r below the window's true rank)
%   exact zeros at the sampled window angles; the optimizer solution provides
%   the well-shaped beam the projection perturbs only slightly.
%   [P7]: candidate null centers are native 2-D (theta, phi) grid points, no
%   1-D cut restriction. Because a 2-D window covers O(width^2) grid points
%   vs O(width) on the old 1-D cut, the projection rank is capped
%   (null_rank_cap) rather than derived purely from the window's numerical
%   rank, to avoid the DOF blowup that destroyed the beam pre-[P7] (P6
%   finding, now handled structurally instead of by restricting to a cut).
%   Arm 1 is the UNCONSTRAINED peak-only solution (null_center_deg = [NaN;NaN]):
%   it serves jammer-off periods and is the peak-gain reference for the P4 gate.
%   Null centers span a (theta, phi) grid at agent_config.null_grid_deg
%   spacing on both axes, excluding points within antijam_config.guard_deg
%   (true spherical separation) of the target. NOTE: grid spacing is O(n^2) in
%   arm count on the full sphere — keep null_grid_deg coarse for 2-D use.
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
%       antijam_config             : struct. Required: theta_s_deg, phi_s_deg,
%                                    guard_deg.
%       agent_config               : struct. Required: null_grid_deg,
%                                    null_width_deg, peak_width_deg,
%                                    null_weight, null_rank_cap.
%       optimizer_config           : struct, per run_optimizer requirements.
%       cache_path                 : char, .mat file for the cached codebook;
%                                    [] disables caching.
%
%   Outputs:
%       codebook : struct with fields
%           W               : (N_el x n_arms) complex weight columns.
%           null_center_deg : (2 x n_arms) arm null centers, row 1 = theta,
%                             row 2 = phi (NaN column for the unconstrained
%                             arm 1).
%           meta            : struct with params (build-parameter echo used
%                             for cache validation) and built (timestamp).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P4, P7].

if nargin < 8
    cache_path = [];
end

theta_s = req_field(antijam_config, 'theta_s_deg', 'antijam');
phi_s   = req_field(antijam_config, 'phi_s_deg',   'antijam');
guard   = req_field(antijam_config, 'guard_deg',   'antijam');
null_grid   = req_field(agent_config, 'null_grid_deg',  'agent');
null_width  = req_field(agent_config, 'null_width_deg', 'agent');
peak_width  = req_field(agent_config, 'peak_width_deg', 'agent');
null_weight = req_field(agent_config, 'null_weight',   'agent');
null_rank_cap = req_field(agent_config, 'null_rank_cap', 'agent');

% ── Candidate null centers on a 2-D (theta, phi) grid, guard cap excluded ──
% Bounded to the actual data grid extent (not the full [0, 180] physical
% range) — many arrays (e.g. a planar array) only have measured pattern data
% over a partial theta span, and a candidate outside that span would yield an
% empty null window below.
theta_cand = min(theta_deg):null_grid:max(theta_deg);
phi_cand   = 0:null_grid:(360 - null_grid / 2);
[TC, PC] = ndgrid(theta_cand, phi_cand);
TC = TC(:).'; PC = PC(:).';
dist = angular_separation_deg(TC, PC, theta_s, phi_s);
keep = dist >= guard;
centers_theta = TC(keep);
centers_phi   = PC(keep);
n_arms = numel(centers_theta) + 1;

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
    'theta_s_deg', theta_s, 'phi_s_deg', phi_s, 'guard_deg', guard, ...
    'null_grid_deg', null_grid, 'null_width_deg', null_width, ...
    'peak_width_deg', peak_width, 'null_weight', null_weight, ...
    'null_rank_cap', null_rank_cap, 'optimizer', optimizer_config, ...
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
peak_dir = struct('type', 'peak', 'theta', theta_s, 'phi', phi_s, ...
                  'width', peak_width, 'weight', 1.0, 'aggregation', 'min');

W = complex(zeros(size(element_patterns_stacked, 1), n_arms));
null_center_deg = NaN(2, n_arms);

fprintf('Codebook: building %d arms (1 unconstrained + %d null sectors)...\n', ...
    n_arms, numel(centers_theta));
result = run_optimizer(element_patterns_stacked, theta_deg, phi_deg, ...
    {peak_dir}, optimizer_config, element_patterns_secondary);
W(:, 1) = result.weights_complex;

for i = 1:numel(centers_theta)
    th_n = centers_theta(i);
    ph_n = centers_phi(i);
    null_dir = struct('type', 'null', 'theta', th_n, 'phi', ph_n, ...
                      'width', null_width, 'weight', null_weight);
    result = run_optimizer(element_patterns_stacked, theta_deg, phi_deg, ...
        {peak_dir, null_dir}, optimizer_config, element_patterns_secondary);
    w_arm = result.weights_complex;

    % Null-space projection polish: exact (rank-capped) zeros on the sampled
    % window angles. [P7]: the full 2-D window is used (no cut restriction),
    % but the projection rank is capped at null_rank_cap to bound the DOF the
    % projection consumes — a 2-D window covers O(width^2) grid points, so
    % using every one of them (as the 1-D cut version safely did) would eat
    % most of the array's DOF and destroy the beam (P6 finding).
    mask = angular_window_mask(theta_deg(:), phi_deg(:), th_n, ph_n, ...
                               null_width, null_width);
    if any(mask(:))
        E_win = window_steering_columns(element_patterns_stacked, mask);
        if ~isempty(element_patterns_secondary)
            E_win = [E_win, window_steering_columns(element_patterns_secondary, mask)]; %#ok<AGROW>
        end
        [U, S, ~] = svd(E_win, 'econ');
        sv = diag(S);
        r_full = sum(sv > 1e-6 * sv(1));
        r = min(r_full, null_rank_cap);
        if r_full > size(W, 1) / 2
            warning('agent_codebook_build:WindowRankHigh', ...
                'Null window at (%g, %g) deg spans rank %d of %d DOF (capped to %d); beam quality may suffer.', ...
                th_n, ph_n, r_full, size(W, 1), r);
        end
        w_arm = w_arm - U(:, 1:r) * (U(:, 1:r)' * w_arm);
    else
        % No grid point falls inside the window (narrow width vs coarse grid
        % resolution) — the optimizer's null directive is the best we can do.
        warning('agent_codebook_build:EmptyNullWindow', ...
            'Null window at (%g, %g) deg contains no grid point; skipping the projection polish for this arm.', ...
            th_n, ph_n);
    end

    W(:, i + 1) = w_arm;
    null_center_deg(:, i + 1) = [th_n; ph_n];
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


function E_win = window_steering_columns(stack, mask)
% Collect the (N_el x M) steering columns of the M grid points where the
% (N_theta x N_phi) window mask is true.
n_el  = size(stack, 1);
flat  = reshape(stack, n_el, []);        % (N_el x N_theta*N_phi), grid flattened
E_win = flat(:, mask(:));
end
