function doa = adapt_music_doa(R_hat, E1, E2, theta_deg, phi_deg, cfg)
% ADAPT_MUSIC_DOA  MUSIC direction-of-arrival estimate of the jammer from R_hat.
%
%   doa = ADAPT_MUSIC_DOA(R_hat, E1, E2, theta_deg, phi_deg, cfg)
%
%   Pure estimator (no state) shared by the predictive Mode C algorithm and the
%   DoA-waterfall plot. Eigendecomposes the running covariance, splits into a
%   fixed 2-source signal subspace (desired + 1 jammer, per locked single-jammer
%   scope) and the remaining noise subspace, and evaluates the MUSIC
%   pseudospectrum
%       P_music(theta,phi) = 1 / (e(theta,phi)' * En * En' * e(theta,phi))
%   over the far-field grid EXCLUDING a guard sector around (theta_s, phi_s).
%   The peak is attributed to the jammer. Jammer PRESENCE is decided from the
%   signal/noise eigengap (2nd eigenvalue vs the noise-floor eigenvalues), NOT
%   from source-order estimation — model order is hardcoded to 2.
%
%   Uses base-MATLAB eig only (R2020a; no Signal Processing / Statistics
%   Toolbox). Steering columns are the measured embedded-element responses from
%   the same E1/E2 stacks the engine uses (mutual coupling included by
%   construction), consistent with sim_engine_step's e_j lookup.
%
%   Inputs:
%       R_hat     : (N_el x N_el) complex Hermitian covariance estimate.
%       E1        : (N_el x N_theta x N_phi) primary-component far-field stack.
%       E2        : (N_el x N_theta x N_phi) secondary component, or [] (1-comp).
%       theta_deg : (1 x N_theta) elevation grid [deg].
%       phi_deg   : (1 x N_phi) azimuth grid [deg].
%       cfg       : struct. Required fields:
%           theta_s_deg, phi_s_deg : desired direction (to build the guard mask).
%           guard_deg              : exclusion radius around theta_s [deg].
%           presence_gap_db        : eigengap threshold to declare "present".
%           doa_stride             : grid subsample factor for the search
%                                    (1 = full grid). Optional; default 1.
%           return_pspec           : logical, fill doa.pspec (waterfall). Opt.
%
%   Outputs:
%       doa : struct with fields
%           present      : logical, jammer detected this step.
%           theta_j_deg  : jammer elevation estimate [deg] (NaN if not present).
%           phi_j_deg    : jammer azimuth estimate [deg]  (NaN if not present).
%           gap_db       : measured signal/noise eigengap [dB] (presence score).
%           sigma_j_sq   : jammer power estimate (linear) from the eigenvalue,
%                          or NaN if not present.
%           pspec        : (N_theta x N_phi) pseudospectrum [dB], or [] unless
%                          cfg.return_pspec is true.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

n_el    = size(R_hat, 1);
n_comp  = 1 + double(~isempty(E2));
n_theta = numel(theta_deg);
n_phi   = numel(phi_deg);
stride  = 1;
if isfield(cfg, 'doa_stride') && ~isempty(cfg.doa_stride)
    stride = max(1, round(cfg.doa_stride));
end
return_pspec = isfield(cfg, 'return_pspec') && ~isempty(cfg.return_pspec) && cfg.return_pspec;

% ── Eigendecomposition of the (Hermitian) covariance ───────────────
% Force exact Hermitian symmetry so eig returns real eigenvalues, then sort
% descending: the top eigenvalues span the signal-plus-jammer subspace.
[V, lam] = eig((R_hat + R_hat') / 2, 'vector');
[lam, ord] = sort(real(lam), 'descend');
V = V(:, ord);

% Source model (single-jammer scope): the desired signal spans n_comp
% eigenvalues and the jammer another n_comp, so the signal subspace has
% dimension 2*n_comp and the noise subspace is everything below it.
n_sig = 2 * n_comp;
if n_sig >= n_el
    error('adapt_music_doa:TooFewElements', ...
        'MUSIC needs N_el > 2*n_comp (= %d); got N_el = %d.', n_sig, n_el);
end
noise_idx   = (n_sig + 1):n_el;
noise_floor = mean(lam(noise_idx));
En = V(:, noise_idx);                          % noise-subspace basis

% ── Presence via the signal/noise eigengap ─────────────────────────
% The first jammer eigenvalue sits at index n_comp+1. When the jammer is off
% it collapses toward the noise floor; a large gap means "jammer present".
gap_db = 10 * log10(lam(n_comp + 1) / noise_floor);
present = gap_db >= cfg.presence_gap_db;

doa = struct('present', present, 'theta_j_deg', NaN, 'phi_j_deg', NaN, ...
    'gap_db', gap_db, 'sigma_j_sq', NaN, 'pspec', []);

% ── MUSIC pseudospectrum over the far-field grid (guard excluded) ──
% (theta, phi) for every grid point, flattened theta-fastest to match the
% engine's idx = (ip-1)*N_theta + it convention (sim_engine_step).
TH = repmat(theta_deg(:),  1, n_phi);
PH = repmat(phi_deg(:).', n_theta, 1);
sep = angular_separation_deg(TH(:), PH(:), cfg.theta_s_deg, cfg.phi_s_deg);
in_guard = sep < cfg.guard_deg;                % exclude the main-beam sector

flat1 = reshape(E1, n_el, []);
flat2 = [];
if n_comp == 2, flat2 = reshape(E2, n_el, []); end

% Candidate grid points: outside the guard, thinned by doa_stride.
cand = find(~in_guard);
cand = cand(1:stride:end);

% Normalized MUSIC metric: ||a||^2 / ||En' a||^2, summed over polarization
% components. Large where a steering column is orthogonal to the noise
% subspace (a source). The ||a||^2 numerator is essential for MEASURED element
% patterns: their column norms vary strongly across the grid (unlike a ULA's
% constant sqrt(N)), so without it the metric spuriously peaks at low-norm
% (e.g. endfire) directions instead of the true source.
A1    = flat1(:, cand);
num   = sum(abs(A1).^2, 1);                    % ||a||^2 per candidate
denom = sum(abs(En' * A1).^2, 1);
if n_comp == 2
    A2    = flat2(:, cand);
    num   = num   + sum(abs(A2).^2, 1);
    denom = denom + sum(abs(En' * A2).^2, 1);
end
pmusic = num ./ denom;

if return_pspec
    pfull = zeros(n_theta * n_phi, 1);         % 0 -> -Inf dB after log
    pfull(cand) = pmusic(:);
    pfull(in_guard) = NaN;                      % blank the excluded sector
    doa.pspec = reshape(10 * log10(pfull), n_theta, n_phi);
end

% ── Peak -> jammer angle (only trusted when the jammer is present) ──
[~, ipk] = max(pmusic);
idx_pk   = cand(ipk);
it = mod(idx_pk - 1, n_theta) + 1;
ip = floor((idx_pk - 1) / n_theta) + 1;
if present
    doa.theta_j_deg = theta_deg(it);
    doa.phi_j_deg   = phi_deg(ip);
    % Jammer power estimate: excess of its eigenvalue over the noise floor,
    % de-embedded by the steering-column norm (per component).
    a_pk = flat1(:, idx_pk);
    if n_comp == 2, a_pk = [a_pk, flat2(:, idx_pk)]; end
    doa.sigma_j_sq = max(lam(n_comp + 1) - noise_floor, 0) / ...
        (sum(vecnorm(a_pk, 2, 1).^2) / n_comp);
end
end
