function state = adapt_spsa_init(spsa_config, w0, seed)
% ADAPT_SPSA_INIT  State for the Mode S SPSA baseline (gradient-free, 2 probes).
%
%   state = ADAPT_SPSA_INIT(spsa_config, w0, seed)
%
%   Warm-startable from any weight vector (uniform weights, a codebook arm, ...).
%   SPSA ASCENDS SINR [dB] on the real parameter vector x = weights_to_x(w)
%   (Milestone-1 re/im interleave; SINR is scale-invariant in w, so the
%   unconstrained parameterization is safe) with gain sequences
%       a_k = a / (A + k)^alpha,   c_k = c / k^gamma
%   (A = Spall's stability constant) and Bernoulli +/-1 simultaneous
%   perturbations delta_k. The ascent step is CLAMPED to norm <= step_max:
%   near a deep null the dB-scale probe difference explodes, and one unclamped
%   step can catapult x_est onto a plateau it never leaves (P3 tuning finding —
%   without the clamp the static-jammer gate fails on most seeds).
%
%   state.w holds the FIRST probe point w(x_est + c_1*delta_1); the runner
%   applies state.w, feeds the observation to adapt_spsa_update, and applies
%   the returned w — each SPSA iteration therefore consumes 2 engine steps.
%
%   Inputs:
%       spsa_config : struct. Required: a, c, alpha, gamma, A, step_max.
%       w0          : (N_el x 1) complex initial weights.
%       seed        : integer, seeds a private RandStream for delta_k.
%
%   Outputs:
%       state : struct with fields
%           x_est     : (2*N_el x 1) current parameter estimate.
%           k         : iteration counter (1-based).
%           phase     : +1 (next obs is the plus-probe) | -1 (minus-probe).
%           delta     : current perturbation direction (+/-1 entries).
%           sinr_plus : stored plus-probe SINR [dB] awaiting the minus probe.
%           cfg       : gain-sequence parameters (a, c, alpha, gamma).
%           w         : (N_el x 1) next weights to apply (a probe point).
%           stream    : private RandStream.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P3].

REQUIRED = {'a', 'c', 'alpha', 'gamma', 'A', 'step_max'};
for i = 1:numel(REQUIRED)
    if ~isfield(spsa_config, REQUIRED{i}) || isempty(spsa_config.(REQUIRED{i}))
        error('adapt_spsa_init:MissingKey', ...
            'Missing required adapt.spsa config key: ''%s''.', REQUIRED{i});
    end
end

x_est  = weights_to_x(w0);
stream = RandStream('mt19937ar', 'Seed', seed);
delta  = 2 * (rand(stream, numel(x_est), 1) > 0.5) - 1;
c_1    = spsa_config.c;   % c_k at k = 1

state = struct( ...
    'x_est',     x_est, ...
    'k',         1, ...
    'phase',     1, ...
    'delta',     delta, ...
    'sinr_plus', NaN, ...
    'cfg',       spsa_config, ...
    'w',         x_to_weights(x_est + c_1 * delta), ...
    'stream',    stream);
end
