function w = adapt_lcmv(R, e_s, diagonal_loading)
% ADAPT_LCMV  MVDR (single-component) / rank-r max-SINR beamformer weights.
%
%   w = ADAPT_LCMV(R, e_s, diagonal_loading)
%
%   Pure function shared by the oracle (analytic R) and the Mode C tracker
%   (estimated R).
%
%   Single component (n_c = 1) — standard MVDR (distortionless toward theta_s):
%       Rl = R + diagonal_loading * eye(N_el)
%       w  = (Rl \ e_s) / (e_s' * (Rl \ e_s)).
%
%   Multiple components (n_c >= 2, i.e. 'total' polarization) — the desired
%   signal is rank-n_c (independent per-component amplitudes), so the correct
%   receiver MAXIMIZES output SINR rather than forcing a distortionless response
%   on every component (which over-constrains the array and collapses the SINR
%   when the components differ greatly in magnitude, e.g. co-pol vs cross-pol).
%   Delegated to adapt_maxsinr; see that function for the derivation. [P8 fix]
%
%   Inputs:
%       R                : (N_el x N_el) complex Hermitian covariance.
%       e_s              : (N_el x n_c) steering column(s) toward theta_s.
%       diagonal_loading : scalar, LINEAR loading level (absolute; callers
%                          convert config diagonal_loading_db relative to the
%                          sigma_n^2 = 1 noise floor). 0 disables loading.
%
%   Outputs:
%       w : (N_el x 1) complex weights.
%
%   Implemented early (during P1): the P1 oracle null-depth test needs it.
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P2, P8].

if size(e_s, 2) == 1
    n_el = size(R, 1);
    Rl   = R + diagonal_loading * eye(n_el);
    A    = Rl \ e_s;
    w    = A / (e_s' * A);                          % MVDR, distortionless
else
    w = adapt_maxsinr(R, e_s, [], diagonal_loading);   % rank-r 'total' pol
end
end
