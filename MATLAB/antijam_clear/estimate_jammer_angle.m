function jammer_estimate = estimate_jammer_angle(covariance, array)
% ══════════════════════════════════════════════════════════════════
% ESTIMATE_JAMMER_ANGLE
% Is a jammer transmitting, and where is it?
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   jammer_estimate = ESTIMATE_JAMMER_ANGLE(covariance, array)
%
%   Two different questions, answered by two different properties of the
%   covariance, and it is worth keeping them straight:
%
%       THE EIGENVALUES SAY WHETHER.  A covariance built from a few strong
%       sources plus noise has a handful of large eigenvalues -- one per source
%       -- sitting above a flat floor of small ones that are noise and nothing
%       else. Counting how many rise above that floor counts the sources: two
%       means signal plus jammer, one means the signal alone. Nothing is
%       hardcoded, so the answer can change during the run, which is what
%       presence detection requires.
%
%       A POWER SCAN SAYS WHERE.  Sweep a unit-length steering vector across
%       the sky and measure the power the array would receive if it pointed
%       there:
%
%           P_beamscan(theta, phi) = e_unit' * R * e_unit
%
%       The largest value outside the guard sector is the jammer. This is the
%       Bartlett (conventional) beamformer response, and it is the estimate that
%       actually places the null.
%
%   MUSIC IS ALSO COMPUTED, AND IS REPORTED RATHER THAN USED.
%
%           P_music(theta, phi) = 1 / || V_noise' * e_unit(theta, phi) ||^2
%
%   MUSIC is the textbook answer and it is sharper on a stationary source,
%   because it measures a NULL in the noise subspace rather than a peak, and a
%   null can be arbitrarily deep. On a MOVING source it is measurably worse
%   here, and the reason is instructive:
%
%       a jammer smeared across ten steps of covariance memory is no longer a
%       point source. Its track occupies a multi-dimensional subspace -- the
%       source count reads 3 on 91 of 120 steps of a drifting run -- so MUSIC
%       makes its noise subspace orthogonal to the WHOLE TRACK, its spectrum
%       goes near-singular along a broad arc, and the peak inside that arc is
%       numerically arbitrary. It sticks, then jumps.
%
%   Measured on spacing0.6, drifting at 0.5 deg/step:
%
%       estimator   median lag   mean |lag|   max |lag|
%       MUSIC          3.50         5.03        13.50     sticks and jumps
%       beamscan       4.50         4.65         6.00     smooth
%
%   The beamscan's lag is simply rate x memory horizon, because a power measure
%   peaks at the power-weighted centre of where the jammer has been -- a
%   predictable error, and one the motion classifier can correct for. MUSIC's
%   is not predictable, and a 13.5 degree excursion cannot be led out.
%
%   On a STATIONARY jammer the two agree exactly, both to 0.00 degrees. MUSIC's
%   real advantage is resolving two closely-spaced sources -- which this
%   single-jammer problem never asks for, while the guard sector already
%   excludes the region near the signal. So the sharper tool is the wrong tool
%   here, and both are reported so that this can be shown rather than claimed.
%
%   TWO THINGS BOTH SCANS MUST GET RIGHT.
%
%   1. UNIT-NORM STEERING VECTORS. Both spectra are only source estimates when
%      e has unit norm. With raw CST vectors the scan multiplies in the array's
%      own gain variation and reports the direction the array hears best rather
%      than the direction the jammer is in. Measured here that error was 77
%      degrees: the jammer was at 55, the scan said 132.
%
%      That is not a tolerable error, it is a fatal one -- nulling at 132 when
%      the jammer is at 55 gives -8.3 dB, WORSE than not adapting at all
%      (-3.6 dB) and 20.4 dB below the correct null. The reason it is called out
%      is that it is SILENT: the spectrum is well formed, the peak is sharp, the
%      answer is confident, and nothing errors. Without ground truth to check
%      against you would ship it. Normalise, and check estimators against truth.
%
%   2. EXCLUDE THE GUARD SECTOR. The wanted signal is one of the strong sources
%      -- it is inside the snapshots -- so both scans find it too, and it is
%      often the stronger peak. Masking the guard sector around the target is
%      exactly what makes "the largest remaining peak" mean "the jammer".
%      Without the mask, Approach 1 confidently and precisely nulls your own
%      signal.
%
%   WHERE THE EIGENVALUE THRESHOLD COMES FROM. With finitely many snapshots the
%   noise eigenvalues do not all sit at the noise power; they spread. For an
%   n-element array and K snapshots the largest reaches about
%
%       (1 + sqrt(n_elements / K))^2   times the true noise power
%
%   (the Marchenko-Pastur upper edge). The threshold is placed there, so it is
%   DERIVED from the array size and the snapshot count rather than tuned, and
%   it is conservative because the covariance also averages over time: it can
%   miss a marginal source, never invent one.
%
%   PRESENCE IS LATE, AND BY MORE THAN IT LOOKS. The threshold sits just above
%   the noise floor, but a strong jammer's eigenvalue starts far above it -- 23.7
%   dB on spacing0.6 at JNR 20 dB -- and the covariance decays it by only
%   -10*log10(lambda) = 0.46 dB per step. It therefore needs about 52 steps to
%   fall below the threshold, so a jammer toggling every 10 steps never appears
%   to switch off at all. The lateness is set by JAMMER STRENGTH, not by the
%   toggle period. This is the accepted cost of leaving out the fast presence
%   covariance -- see README.md, "Deliberate omissions".
%
%   Inputs:
%       covariance : (n_elements x n_elements) Hermitian, from
%                    sample_covariance. Units: power.
%       array      : struct from make_array. Supplies the element patterns, the
%                    angle grids, K = snapshots_per_step (which sets the
%                    eigenvalue threshold) and the profile, which supplies the
%                    target direction and guard_deg.
%
%   Outputs:
%       jammer_estimate : struct with fields
%           is_feasible      : logical, false when the array is too small to
%                              tell how many sources it is hearing.
%           is_present       : logical, a second source was found.
%           theta_deg        : jammer elevation from the BEAMSCAN -- the
%                              estimate that places the null. NaN if absent.
%           phi_deg          : jammer azimuth from the beamscan. Units: degrees.
%           music_theta_deg  : the same from MUSIC, for comparison only.
%           music_phi_deg    : ditto. Units: degrees.
%           n_sources        : eigenvalues above the noise floor.
%           peak_ratio_db    : beamscan peak over the median of the searched
%                              region -- a sharpness figure. Units: dB.
%           reason           : explanation when infeasible or absent.

% ────────────────────────── CONSTANTS ─────────────────────────────

% Guards a division and a logarithm; far below any physical value here.
TINY = 1e-30;

% A source must exceed the noise floor by at least the Marchenko-Pastur edge
% (computed below) times this margin. 1.0 means exactly at the predicted edge --
% no fudge factor. It is named so it is visible that there is not one.
EIGENVALUE_MARGIN = 1.0;

% Signal plus jammer. Fewer than this means the jammer is not transmitting.
SOURCES_WHEN_JAMMED = 2;

element_patterns = array.element_patterns;
theta_deg        = array.theta_deg(:);
phi_deg          = array.phi_deg(:);
profile          = array.profile;
n_snapshots      = array.snapshots_per_step;
n_elements       = size(covariance, 1);

jammer_estimate = struct( ...
    'is_feasible',     false, ...
    'is_present',      false, ...
    'theta_deg',       NaN, ...
    'phi_deg',         NaN, ...
    'music_theta_deg', NaN, ...
    'music_phi_deg',   NaN, ...
    'n_sources',       NaN, ...
    'peak_ratio_db',   NaN, ...
    'reason',          '');

% ────────────────────────── FEASIBILITY ───────────────────────────

% At least one eigenvector must be left over after the sources, or there is no
% noise floor to compare anything against and presence cannot be decided. With
% a signal and a jammer that means three elements. A two-element array can form
% a null but cannot reliably find the jammer, so it correctly declines -- and
% saying so is the right output, not an exception or a fabricated angle.
if n_elements < SOURCES_WHEN_JAMMED + 1
    jammer_estimate.reason = sprintf( ...
        ['This array has %d element(s). Telling a jammer from the noise needs ' ...
         'more elements than sources, i.e. at least %d for a signal and a ' ...
         'jammer. Direction finding is not possible here.'], ...
        n_elements, SOURCES_WHEN_JAMMED + 1);
    return
end
jammer_estimate.is_feasible = true;

% ────────────────────────── HOW MANY SOURCES? ─────────────────────

[eigenvectors, eigenvalues] = eig(covariance, 'vector');
eigenvalues = real(eigenvalues);

[eigenvalues, order] = sort(eigenvalues, 'descend');
eigenvectors         = eigenvectors(:, order);

% The same noise floor the diagonal loading is scaled to -- one definition,
% two callers, no way for them to disagree.
noise_floor      = noise_floor_power(covariance);
sampling_edge    = (1 + sqrt(n_elements / n_snapshots)) ^ 2;
source_threshold = noise_floor * sampling_edge * EIGENVALUE_MARGIN;

n_sources = sum(eigenvalues > source_threshold);

% At least one eigenvector must remain for the noise subspace.
n_sources = max(1, min(n_sources, n_elements - 1));
jammer_estimate.n_sources = n_sources;

if n_sources < SOURCES_WHEN_JAMMED
    jammer_estimate.reason = sprintf( ...
        ['Only %d source above the noise floor: the wanted signal alone. ' ...
         'The jammer is not transmitting (or is too weak to resolve).'], n_sources);
    return
end
jammer_estimate.is_present = true;

% ────────────────────────── THE SEARCH GRID ───────────────────────

n_theta  = numel(theta_deg);
n_phi    = numel(phi_deg);
manifold = reshape(element_patterns, n_elements, n_theta * n_phi);

% Unit norm per direction -- see note 1 in the header. Shared by both scans.
manifold_unit = manifold ./ max(vecnorm(manifold, 2, 1), TINY);

% Only directions far enough from the target to be a jammer at all.
separation_deg = angular_separation_grid(theta_deg, phi_deg, ...
                                         profile.signal_theta_deg, profile.signal_phi_deg);
is_searchable  = separation_deg > profile.guard_deg;

if ~any(is_searchable(:))
    jammer_estimate.is_present = false;
    jammer_estimate.reason = sprintf( ...
        ['The guard sector of %.1f deg covers the entire sphere for this ' ...
         '(array, target) pair, so there is nowhere a jammer could be nulled. ' ...
         'The beam is too broad for this method to apply here.'], profile.guard_deg);
    return
end

% ────────────────────────── SCAN 1: BEAMSCAN (USED) ───────────────

% Power received if the array pointed each way: e_unit' * R * e_unit.
beamscan_spectrum = real(sum(conj(manifold_unit) .* (covariance * manifold_unit), 1));
beamscan_spectrum = reshape(beamscan_spectrum, n_theta, n_phi);

[jammer_estimate.theta_deg, jammer_estimate.phi_deg, beamscan_peak] = ...
    masked_peak(beamscan_spectrum, is_searchable, theta_deg, phi_deg);

% Sharpness against the typical level of the searched region. A genuine source
% stands well above it; a flat spectrum means nothing in particular was found.
background = median(beamscan_spectrum(is_searchable));
jammer_estimate.peak_ratio_db = 10 * log10(max(beamscan_peak, TINY) / max(background, TINY));

% ────────────────────────── SCAN 2: MUSIC (REPORTED) ──────────────

noise_subspace   = eigenvectors(:, n_sources + 1 : end);
noise_projection = noise_subspace' * manifold_unit;
noise_fraction   = sum(abs(noise_projection) .^ 2, 1);

music_spectrum = reshape(1 ./ max(noise_fraction, TINY), n_theta, n_phi);

[jammer_estimate.music_theta_deg, jammer_estimate.music_phi_deg] = ...
    masked_peak(music_spectrum, is_searchable, theta_deg, phi_deg);

jammer_estimate.reason = sprintf( ...
    'Beamscan peak at (%.1f, %.1f) deg, %.1f dB above background; %d sources.', ...
    jammer_estimate.theta_deg, jammer_estimate.phi_deg, ...
    jammer_estimate.peak_ratio_db, n_sources);
end


% ────────────────────────── HELPERS ───────────────────────────────

function [peak_theta_deg, peak_phi_deg, peak_value] = masked_peak(spectrum, is_searchable, ...
                                                                  theta_deg, phi_deg)
% Largest value of a spectrum within the searchable region, and where it is.
spectrum(~is_searchable) = -Inf;

[peak_value, flat_index] = max(spectrum(:));
[theta_index, phi_index] = ind2sub(size(spectrum), flat_index);

peak_theta_deg = theta_deg(theta_index);
peak_phi_deg   = phi_deg(phi_index);
end


function separation_deg = angular_separation_grid(theta_deg, phi_deg, ...
                                                  reference_theta_deg, reference_phi_deg)
% Great-circle angle from a reference direction to every point of the grid.
%
% Using the true angle on the sphere, rather than a difference in theta alone,
% matters near the poles, where a large change in phi is a small change in
% direction.
%
%   cos(separation) = cos(t) cos(t0) + sin(t) sin(t0) cos(phi - phi0)

theta_rad           = deg2rad(theta_deg(:));         % (n_theta x 1)
phi_rad             = deg2rad(phi_deg(:)).';         % (1 x n_phi)
reference_theta_rad = deg2rad(reference_theta_deg);
reference_phi_rad   = deg2rad(reference_phi_deg);

% [MATLAB] implicit expansion gives the (n_theta x n_phi) grid directly.
cos_separation = cos(theta_rad) * cos(reference_theta_rad) ...
                 + sin(theta_rad) * sin(reference_theta_rad) ...
                   .* cos(phi_rad - reference_phi_rad);

separation_deg = rad2deg(acos(min(max(cos_separation, -1.0), 1.0)));
end
