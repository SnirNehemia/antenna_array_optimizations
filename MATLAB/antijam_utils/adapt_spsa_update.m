function [w, state] = adapt_spsa_update(state, obs)
% ADAPT_SPSA_UPDATE  Mode S update: alternate SPSA probes, step on probe pairs.
%
%   [w, state] = ADAPT_SPSA_UPDATE(state, obs)
%
%   Consumes obs.sinr_db ONLY (Mode S contract). Each call returns the next
%   probe point; every second call additionally applies the SPSA ascent step:
%       phase +1: store obs.sinr_db as the pending pair's plus-probe value,
%                 return w(x_est - c_k * delta)             (minus probe next)
%       phase -1: g_hat = (sinr_plus - obs.sinr_db) / (2 * c_k) * delta
%                 (delta has +/-1 entries, so ./delta == .*delta),
%                 x_est <- x_est + a_k * g_hat with the step norm clamped to
%                 cfg.step_max, k <- k + 1, draw new delta,
%                 return w(x_est + c_k * delta)              (plus probe next)
%   Gains: a_k = a / (A + k)^alpha, c_k = c / k^gamma (see adapt_spsa_init on
%   why the clamp and stability constant are load-bearing).
%
%   Inputs:
%       state : struct from adapt_spsa_init.
%       obs   : observation struct from sim_engine_step (Mode S).
%
%   Outputs:
%       w     : (N_el x 1) complex weights to apply next step (a probe point).
%       state : updated state.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P3].

cfg = state.cfg;
c_k = cfg.c / state.k^cfg.gamma;

if state.phase == 1
    % Plus probe observed; issue the minus probe of the same pair.
    state.sinr_plus = obs.sinr_db;
    state.phase     = -1;
    w = x_to_weights(state.x_est - c_k * state.delta);
else
    % Minus probe observed; take the clamped ascent step, open the next pair.
    a_k   = cfg.a / (cfg.A + state.k)^cfg.alpha;
    g_hat = (state.sinr_plus - obs.sinr_db) / (2 * c_k) * state.delta;
    step  = a_k * g_hat;
    step_norm = norm(step);
    if step_norm > cfg.step_max
        step = step * (cfg.step_max / step_norm);
    end
    state.x_est = state.x_est + step;

    state.k     = state.k + 1;
    state.delta = 2 * (rand(state.stream, numel(state.x_est), 1) > 0.5) - 1;
    state.phase = 1;
    c_k = cfg.c / state.k^cfg.gamma;
    w = x_to_weights(state.x_est + c_k * state.delta);
end
state.w = w;
end
