function [w, state] = adapt_tracking_update(state, obs)
% ADAPT_TRACKING_UPDATE  Mode C update: fold in snapshots, recompute LCMV weights.
%
%   [w, state] = ADAPT_TRACKING_UPDATE(state, obs)
%
%   Consumes obs.snapshots ONLY (Mode C contract; obs.sinr_db is ignored so the
%   tracker never depends on the scalar channel). For each snapshot column x:
%       R_hat <- lambda * R_hat + (1 - lambda) * (x * x')
%   then w = adapt_lcmv(R_hat, e_s, loading). Note the snapshots contain the
%   desired signal (engine snapshot model), so this is an MPDR-style tracker —
%   diagonal loading is the guard against finite-sample signal self-nulling.
%
%   Inputs:
%       state : struct from adapt_tracking_init.
%       obs   : observation struct from sim_engine_step (Mode C).
%
%   Outputs:
%       w     : (N_el x 1) complex weights to apply next step.
%       state : updated state (R_hat, w).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P2].

if isempty(obs.snapshots)
    error('adapt_tracking_update:NoSnapshots', ...
        'obs.snapshots is empty — the covariance tracker requires Mode C.');
end

X = obs.snapshots;
for i = 1:size(X, 2)
    x = X(:, i);
    state.R_hat = state.lambda * state.R_hat + (1 - state.lambda) * (x * x');
end

w = adapt_lcmv(state.R_hat, state.e_s, state.loading);
state.w = w;
end
