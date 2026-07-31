function w = adapt_lcmv_null(R, e_s, e_null, diagonal_loading)
% ADAPT_LCMV_NULL  Multi-constraint LCMV: distortionless at theta_s, null at theta_j.
%
%   w = ADAPT_LCMV_NULL(R, e_s, e_null, diagonal_loading)
%
%   Generalizes adapt_lcmv (which stays unmodified as the P2 single-constraint
%   MVDR) to a hard null constraint. With constraint matrix C = [e_s, e_null]
%   and response vector g = [ones(n_s,1); zeros(n_null,1)]:
%       Rl = R + diagonal_loading * eye(N_el)
%       w  = Rl \ C * ((C' * (Rl \ C)) \ g)
%   The null is a DETERMINISTIC constraint on w, independent of R — so it
%   pre-nulls a predicted angle even when that angle currently carries no energy
%   (the anticipatory / pre-turn-on behavior P8 relies on). When e_null is empty
%   this reduces exactly to adapt_lcmv(R, e_s, diagonal_loading) (asserted in the
%   P8 unit test).
%
%   Inputs:
%       R                : (N_el x N_el) complex Hermitian covariance. Use the
%                          estimated R_hat when the jammer is present, or eye
%                          (quiescent) when pre-nulling a currently-silent angle.
%       e_s              : (N_el x n_s) steering column(s) toward theta_s
%                          (per-component; distortionless response).
%       e_null           : (N_el x n_null) steering column(s) toward the
%                          (predicted) jammer angle, or [] for no null.
%       diagonal_loading : scalar, LINEAR loading level (absolute; callers
%                          convert diagonal_loading_db relative to the
%                          sigma_n^2 = 1 noise floor). 0 disables loading.
%
%   Outputs:
%       w : (N_el x 1) complex weights.
%
%   For 'total' polarization (n_c >= 2 signal columns) a distortionless
%   constraint on every component is over-constrained; the desired signal is
%   rank-n_c, so max-SINR-with-hard-null (adapt_maxsinr) is used instead. [P8 fix]
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

if size(e_s, 2) >= 2
    w = adapt_maxsinr(R, e_s, e_null, diagonal_loading);   % rank-r 'total' pol
    return;
end

n_el = size(R, 1);

% Constraint matrix C and response vector g: distortionless (1) on the single
% e_s column, null (0) on every e_null column. e_null = [] -> plain MVDR.
C = [e_s, e_null];
g = [1; zeros(size(e_null, 2), 1)];

% Loaded-covariance LCMV solution (same form as adapt_lcmv, extended to g):
%   w = Rl^-1 C (C' Rl^-1 C)^-1 g.
% pinv on the small (n_constraints x n_constraints) gram is used instead of \:
% it is identical to \ when the constraints are independent, but degrades
% gracefully (min-norm, no warning) if a predicted null direction happens to
% coincide with the steering direction (e.g. a ULA's theta/(180-theta)
% endfire ambiguity) and makes C rank-deficient.
Rl = R + diagonal_loading * eye(n_el);
A  = Rl \ C;                                   % (N_el x n_constraints)
w  = A * (pinv(C' * A) * g);
end
