function w = adapt_lcmv(R, e_s, diagonal_loading)
% ADAPT_LCMV  LCMV/MVDR weights from a covariance and steering constraint(s).
%
%   w = ADAPT_LCMV(R, e_s, diagonal_loading)
%
%   Pure function shared by the oracle (analytic R) and the Mode C tracker
%   (estimated R). With C = e_s (N_el x n_c constraint columns) and
%   f = ones(n_c, 1) (distortionless response per component):
%       Rl = R + diagonal_loading * eye(N_el)
%       w  = Rl \ C * ((C' * (Rl \ C)) \ f)
%   For a single component (n_c = 1) this reduces to MVDR:
%       w = Rl \ e_s / (e_s' * (Rl \ e_s)).
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
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P2].

n_el = size(R, 1);
Rl   = R + diagonal_loading * eye(n_el);
A    = Rl \ e_s;                                   % (N_el x n_c)
w    = A * ((e_s' * A) \ ones(size(e_s, 2), 1));
end
