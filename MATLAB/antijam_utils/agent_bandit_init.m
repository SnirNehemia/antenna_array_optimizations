function state = agent_bandit_init(agent_config, codebook, seed)
% AGENT_BANDIT_INIT  State for the non-stationary bandit over codebook arms.
%
%   state = AGENT_BANDIT_INIT(agent_config, codebook, seed)
%
%   Primary method: DISCOUNTED THOMPSON SAMPLING with a Gaussian reward model
%   per arm. Sufficient statistics (discounted count n_i and reward sum s_i)
%   decay by agent_config.discount on every update, so unplayed arms lose
%   evidence and get re-explored — the mechanism that tracks a moving jammer.
%   Posterior sample per arm: N(s_i/n_i, SIGMA_TILDE^2 / n_i). Arms start with
%   an OPTIMISTIC prior (mean = CLIP_HI, tiny count), so the first ~n_arms
%   updates sweep every arm once before exploitation begins.
%
%   Alternative (agent_config.method = 'swucb'): sliding-window UCB over the
%   last agent_config.window plays: mean of in-window rewards per arm
%   + SIGMA_TILDE * sqrt(2*log(t_w)/n_i), unplayed-in-window arms first.
%
%   Reward = obs.sinr_db clipped to [CLIP_LO, CLIP_HI] = [-10, +40] dB
%   (plan Section 7 decision #3, locked in P0).
%
%   Inputs:
%       agent_config : struct. Required: method ('thompson'|'swucb'), discount,
%                      window.
%       codebook     : struct from agent_codebook_build.
%       seed         : integer, seeds a private RandStream (posterior sampling).
%
%   Outputs:
%       state : struct with fields
%           method    : 'thompson' | 'swucb'.
%           W         : (N_el x n_arms) arm weights (copied from codebook).
%           arm_index : currently applied arm (1-based; starts at 1, the
%                       unconstrained arm).
%           w         : (N_el x 1) weights of the current arm.
%           n, s      : (1 x n_arms) discounted count / reward sum (thompson).
%           win_arm, win_reward, win_pos : circular play history (swucb).
%           cfg       : discount / window / clip bounds / sigma_tilde.
%           stream    : private RandStream.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P5].

REQUIRED = {'method', 'discount', 'window', 'sigma_tilde_db'};
for i = 1:numel(REQUIRED)
    if ~isfield(agent_config, REQUIRED{i}) || isempty(agent_config.(REQUIRED{i}))
        error('agent_bandit_init:MissingKey', ...
            'Missing required agent config key: ''%s''.', REQUIRED{i});
    end
end
if ~any(strcmp(agent_config.method, {'thompson', 'swucb'}))
    error('agent_bandit_init:BadMethod', ...
        'Unknown bandit method ''%s'' (expected ''thompson'' or ''swucb'').', ...
        agent_config.method);
end

CLIP_LO = -10;      % dB — plan decision #3 (locked P0)
CLIP_HI =  40;
% Posterior/confidence scale [dB] (config agent.sigma_tilde_db). Small: the
% pattern-level reward is near-deterministic, so discrimination is sharp
% after one pull; re-exploration under drift comes from the DISCOUNT decaying
% n_i. It also sets the exploration-overhead rate: far arms get resampled
% when sigma_tilde/sqrt(n_i) grows to their reward gap. (P5 sweep: 5 dB ->
% 40% static identification, 1 dB -> ~95%; P6 real-data sweep drove it lower
% still to cut steady-state exploration.)
SIGMA_TILDE = agent_config.sigma_tilde_db;
PRIOR_COUNT = 1e-2; % optimistic prior weight

n_arms = size(codebook.W, 2);
w_len  = agent_config.window;

state = struct();
state.method    = agent_config.method;
state.W         = codebook.W;
state.arm_index = 1;
state.w         = codebook.W(:, 1);
state.n         = PRIOR_COUNT * ones(1, n_arms);
state.s         = PRIOR_COUNT * CLIP_HI * ones(1, n_arms);   % optimistic mean
state.win_arm    = zeros(1, w_len);    % 0 = empty slot
state.win_reward = zeros(1, w_len);
state.win_pos    = 0;
state.cfg = struct('discount', agent_config.discount, 'window', w_len, ...
                   'clip_lo', CLIP_LO, 'clip_hi', CLIP_HI, ...
                   'sigma_tilde', SIGMA_TILDE);
state.stream = RandStream('mt19937ar', 'Seed', seed);
end
