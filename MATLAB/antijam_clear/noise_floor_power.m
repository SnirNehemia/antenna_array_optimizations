function noise_power = noise_floor_power(covariance)
% ══════════════════════════════════════════════════════════════════
% NOISE_FLOOR_POWER
% How much of what the array hears is just receiver noise.
%
% Part of: Antenna Array Pattern Optimization Tool -- anti-jam (clear rewrite)
% ══════════════════════════════════════════════════════════════════
%
%   noise_power = NOISE_FLOOR_POWER(covariance)
%
%   A covariance built from a few strong sources plus noise has a characteristic
%   eigenvalue profile: a handful of large values, one per source, and then a
%   flat floor of small ones that are noise and nothing else. This returns the
%   level of that floor.
%
%   It is read from the SMALLER HALF of the eigenvalues, which are noise by
%   construction whenever there are only a couple of sources, and the MEDIAN is
%   used rather than the mean so that one source leaking into the set cannot
%   drag the estimate up.
%
%   Two callers need it and they need the same answer:
%       estimate_jammer_angle -- to decide which eigenvalues are sources;
%       max_sinr_weights      -- to scale the diagonal loading.
%
%   WHY THE NOISE FLOOR RATHER THAN THE TRACE. The obvious alternative is to
%   scale things by trace(R)/n, the average power. That is wrong for loading:
%   the trace is dominated by the jammer, so a stronger jammer would raise your
%   own diagonal loading and blunt the very null you need against it. The noise
%   floor is the one part of the covariance that does not move when the jammer
%   does, which is exactly the property a reference level should have.
%
%   Inputs:
%       covariance : (n_elements x n_elements) Hermitian. Units: power.
%
%   Outputs:
%       noise_power : scalar estimate of the per-element noise power.
%                     Units: power (same as the covariance).

n_elements  = size(covariance, 1);
eigenvalues = sort(real(eig(covariance)), 'descend');

% Drop the upper half, which may contain sources.
lower_half = eigenvalues(ceil(n_elements / 2) + 1 : end);
if isempty(lower_half)
    % One- or two-element arrays have no "half" to speak of. The smallest
    % eigenvalue is still the best available estimate of the floor.
    lower_half = eigenvalues(end);
end

noise_power = median(lower_half);

% A covariance always has noise in it -- this cannot legitimately be zero, and a
% zero here would make the diagonal loading vanish silently.
if ~isfinite(noise_power) || noise_power <= 0
    error('noise_floor_power:NoNoiseFloor', ...
        ['Estimated noise floor is %g, which is not a usable power. The ' ...
         'covariance is singular or non-finite; check the snapshot block.'], ...
        noise_power);
end
end
