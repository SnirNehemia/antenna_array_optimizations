function w = adapt_maxsinr(R, e_s, e_null, diagonal_loading)
% ADAPT_MAXSINR  Rank-r max-output-SINR beamformer, with an optional hard null.
%
%   w = ADAPT_MAXSINR(R, e_s, e_null, diagonal_loading)
%
%   Maximizes the output SINR of a rank-r desired signal:
%       maximize  w' R_s w / w' Rl w      subject to   w' e_null = 0,
%   with  R_s = e_s e_s'  (summed over the r desired-signal steering columns)
%   and   Rl  = R + diagonal_loading * I.
%
%   r = size(e_s, 2): r = 1 for a single polarization component, r = 2 for
%   'total' polarization (co-pol + cross-pol, independent per-component
%   amplitudes). This is the CORRECT max-SINR receiver for a rank-r signal and
%   replaces imposing a distortionless constraint on EVERY component, which is
%   over-constrained: forcing w'e_copol = w'e_cross = 1 ties two very different
%   spatial signatures to the same gain and, when the components differ greatly
%   in magnitude, wastes array gain and collapses the SINR. For r = 1 and
%   e_null = [] this reduces to MVDR (w propto Rl^-1 e_s).
%
%   Closed form (no null): the optimum lies in span(Rl^-1 e_s). With
%   A = Rl^-1 e_s (N x r) and the r x r Hermitian Gram G = e_s' A, the SINR is
%   maximized by w = A v where v is the principal eigenvector of G. With a hard
%   null, the same is solved inside an orthonormal basis Q of {w : e_null' w = 0}
%   (so w = Q z, w' e_null = 0 exactly). SINR and the radiation pattern are
%   scale-invariant, so w is returned unit-norm.
%
%   Inputs:
%       R                : (N_el x N_el) complex Hermitian covariance.
%       e_s              : (N_el x r) desired-signal steering column(s).
%       e_null           : (N_el x r_null) steering column(s) to null, or [].
%       diagonal_loading : scalar, LINEAR loading added to R. 0 disables it.
%
%   Outputs:
%       w : (N_el x 1) complex weights (unit-norm).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P8].

n_el = size(R, 1);
Rl   = R + diagonal_loading * eye(n_el);

if isempty(e_null)
    A = Rl \ e_s;                          % (N x r) = Rl^-1 e_s
    G = e_s' * A;                          % (r x r) Hermitian Gram
    v = principal_eigvec(G);
    w = A * v;
else
    Q  = null(e_null');                    % (N x N-r_null): {w : e_null' w = 0}
    Es = Q' * e_s;                         % steering in the null-constrained subspace
    A  = (Q' * Rl * Q) \ Es;
    G  = Es' * A;
    z  = A * principal_eigvec(G);
    w  = Q * z;
end

w = w / norm(w);
end


function v = principal_eigvec(G)
% Eigenvector of the Hermitian matrix G for its largest eigenvalue.
[V, D] = eig((G + G') / 2);
[~, ix] = max(real(diag(D)));
v = V(:, ix);
end
