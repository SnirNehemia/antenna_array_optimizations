function [w, state] = agent_bandit_update(state, obs)
% AGENT_BANDIT_UPDATE  Mode S update: credit reward to current arm, pick next.
%
%   [w, state] = AGENT_BANDIT_UPDATE(state, obs)
%
%   Consumes obs.sinr_db ONLY (Mode S contract). Credits the clipped reward to
%   state.arm_index, then selects the next arm:
%       thompson : discount ALL arms' statistics, credit the played arm, draw
%                  theta_i ~ N(s_i/n_i, sigma_tilde^2/n_i) per arm, play the
%                  argmax. Discounting decays unplayed arms' counts, inflating
%                  their sampling variance -> periodic re-exploration.
%       swucb    : record the play in the circular window, play the argmax of
%                  in-window mean + sigma_tilde*sqrt(2*log(t_w)/n_i); arms
%                  with no in-window plays come first.
%   Returns w = state.W(:, next_arm).
%
%   Inputs:
%       state : struct from agent_bandit_init.
%       obs   : observation struct from sim_engine_step (Mode S).
%
%   Outputs:
%       w     : (N_el x 1) complex weights (the selected arm's column).
%       state : updated state (stats, arm_index, w).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P5].

cfg = state.cfg;
r   = min(max(obs.sinr_db, cfg.clip_lo), cfg.clip_hi);
j   = state.arm_index;

if strcmp(state.method, 'thompson')
    % Discount everyone, credit the played arm.
    state.n = cfg.discount * state.n;
    state.s = cfg.discount * state.s;
    state.n(j) = state.n(j) + 1;
    state.s(j) = state.s(j) + r;
    % Gaussian posterior draw per arm.
    mu      = state.s ./ state.n;
    sig     = cfg.sigma_tilde ./ sqrt(state.n);
    samples = mu + sig .* randn(state.stream, 1, numel(mu));
    [~, next_arm] = max(samples);
else % swucb
    state.win_pos = mod(state.win_pos, cfg.window) + 1;
    state.win_arm(state.win_pos)    = j;
    state.win_reward(state.win_pos) = r;
    n_arms = size(state.W, 2);
    ucb = zeros(1, n_arms);
    t_w = nnz(state.win_arm);
    for i = 1:n_arms
        in_win = state.win_arm == i;
        n_i = nnz(in_win);
        if n_i == 0
            ucb(i) = Inf;                          % unplayed in window: force
        else
            ucb(i) = mean(state.win_reward(in_win)) ...
                   + cfg.sigma_tilde * sqrt(2 * log(t_w) / n_i);
        end
    end
    [~, next_arm] = max(ucb);
end

state.arm_index = next_arm;
w = state.W(:, next_arm);
state.w = w;
end
