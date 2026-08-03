function [obs, sim_state] = sim_engine_step(sim_state, w, f_notch_hz)
% SIM_ENGINE_STEP  One closed-loop step: apply weights, observe SINR (+snapshots).
%
%   [obs, sim_state] = SIM_ENGINE_STEP(sim_state, w)
%   [obs, sim_state] = SIM_ENGINE_STEP(sim_state, w, f_notch_hz)
%
%   THE ONLY observation channel available to adaptation algorithms. Advances
%   the scenario one step (k -> k+1) and evaluates the applied weights against
%   the true jammer state, which stays hidden inside sim_state.
%
%   SINR (pattern-level, per Section 2 of antijam_milestone_plan.md):
%       P_sig = sigma_s^2 * sum_c |w' * e_c(theta_s)|^2 / n_comp
%       P_jam = sigma_j^2 * sum_c |w' * e_c(theta_j)|^2 / n_comp   (0 when off)
%       SINR  = P_sig / (P_jam + sigma_n^2 * norm(w)^2)
%   where c runs over the 1 or 2 polarization components (equal power split).
%
%   Snapshot model (Mode C):
%       x(t) = sum_c e_c(theta_s) s_c(t) + sum_c e_c(theta_j) j_c(t) + n(t)
%   with s_c, n i.i.d. circular complex Gaussian drawn from sim_state.stream as
%   sqrt(sigma^2/2) * (randn + 1i*randn) — temporally white, i.e. a wideband
%   desired signal and a white noise floor across the captured band.
%
%   [P10] WAVEFORM LAYER (opt-in, sim_state.use_waveform):
%   the jammer amplitudes j_c(t) become a CW TONE at the scenario's carrier
%   offset (sim_jammer_waveform) instead of white Gaussian, giving the snapshot
%   block a spectral line for adapt_freq_estimate to find. This leaves the
%   SPATIAL covariance E[x x'] unchanged, so every Mode C algorithm — which sees
%   only (X*X')/K — is unaffected. With use_waveform false the generator below
%   is byte-identical to the pre-[P10] path under the same seed.
%
%   [P10] RF NOTCH (opt-in, sim_state.notch): a linear time-invariant notch is
%   identical on every element and therefore COMMUTES with the beamformer, so it
%   enters the pattern-level SINR as exactly two scalars from sim_notch_response:
%       P_sig   *= (1 - loss_frac)      % wideband signal loses the notched slice
%       P_noise *= (1 - loss_frac)      % white noise loses the same slice
%       P_jam   *= |H(f_j_true)|^2      % notch centred on the ESTIMATE, evaluated
%                                       % at TRUTH -> estimation error costs
%                                       % rejection, which is the point
%   The notch is commanded by the algorithm through f_notch_hz, exactly as the
%   spatial response is commanded through w. notch.mode 'oracle' ignores the
%   command and uses truth (upper bound); 'off' disables the notch; 'adaptive'
%   honours the command and is the real system. An uncommanded notch (f_notch_hz
%   omitted or NaN) is simply not engaged and costs nothing.
%
%   Inputs:
%       sim_state  : struct from sim_engine_init (advanced state returned).
%       w          : (N_el x 1) complex weights to apply this step.
%       f_notch_hz : [P10] optional commanded notch centre, as an offset from
%                    the array centre frequency [Hz]. NaN / omitted -> notch not
%                    engaged this step.
%
%   Outputs:
%       obs : struct with fields
%           sinr_db   : scalar, output SINR of w at the current true state [dB],
%                       including the notch's effect when one is engaged.
%           snapshots : (N_el x K) complex, K = sim_state.n_snapshots — Mode C
%                       only; [] in Mode S. Mode-S algorithms MUST NOT touch it.
%           fs_hz     : [P10] snapshot sample rate [Hz] accompanying snapshots
%                       (Mode C with the waveform layer on); [] otherwise.
%       sim_state : updated state (k advanced; stream advances in place — it is
%                   a handle object). sim_state.last carries truth-side notch
%                   diagnostics (sinr_no_notch_db, f_notch_hz, notch_rejection_db,
%                   loss_frac) for closed_loop_run / kpi_evaluate — these are
%                   NOT part of obs and must never reach an algorithm.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P1, P10].

if nargin < 3
    f_notch_hz = NaN;
end

k = sim_state.k + 1;
scn = sim_state.scenario;
if k > numel(scn.t_s)
    error('sim_engine_step:EndOfScenario', ...
        'Scenario ''%s'' has %d steps; step %d requested.', ...
        scn.id, numel(scn.t_s), k);
end
sim_state.k = k;

w = w(:);
n_el   = size(sim_state.E1, 1);
n_comp = 1 + double(~isempty(sim_state.E2));

% True jammer steering column(s) at the current step.
[it_j, ip_j] = nearest_index_2d(sim_state.theta_deg, sim_state.phi_deg, ...
    scn.theta_j_deg(k), scn.phi_j_deg(k));
idx_j = (ip_j - 1) * numel(sim_state.theta_deg) + it_j;
e_j   = sim_state.E1(:, idx_j);
if n_comp == 2
    e_j = [e_j, sim_state.E2(:, idx_j)];
end
if scn.jammer_on(k)
    sigma_j_sq = 10^(scn.jn_ratio_db(k) / 10);
else
    sigma_j_sq = 0;
end

% ── [P10] Notch state for this step ───────────────────────────────
% h2_jam scales the jammer power; loss_frac is the wideband insertion loss on
% the desired signal and the noise. Both are 1 / 0 when no notch is engaged.
f_j_true   = 0;
if sim_state.use_waveform
    f_j_true = scn.f_j_hz(k);
end
[h2_jam, loss_frac, f_notch_used] = notch_state(sim_state, f_notch_hz, f_j_true);

% ── Pattern-level SINR ────────────────────────────────────────────
p_sig   = sim_state.sigma_s_sq * sum(abs(w' * sim_state.e_s).^2) / n_comp;
p_jam   = sigma_j_sq           * sum(abs(w' * e_j).^2)           / n_comp;
p_noise = sim_state.sigma_n_sq * real(w' * w);

sinr_no_notch = p_sig / (p_jam + p_noise);
sinr_notched  = (p_sig * (1 - loss_frac)) / ...
                (p_jam * h2_jam + p_noise * (1 - loss_frac));

obs = struct('sinr_db', 10 * log10(sinr_notched), 'snapshots', [], 'fs_hz', []);

% Truth-side diagnostics: reporting only, never observable by an algorithm.
sim_state.last = struct( ...
    'sinr_no_notch_db',   10 * log10(sinr_no_notch), ...
    'f_notch_hz',         f_notch_used, ...
    'notch_rejection_db', -10 * log10(h2_jam), ...
    'loss_frac',          loss_frac);

% ── Mode C snapshots ──────────────────────────────────────────────
if strcmp(sim_state.mode, 'C')
    st = sim_state.stream;   % handle object: draws advance it in place
    K  = sim_state.n_snapshots;
    % Accumulation order (signal -> jammer -> noise) is load-bearing: it keeps
    % the pre-[P10] path bit-for-bit reproducible, floating-point rounding
    % included. Do not reorder.
    x   = sim_state.e_s * cgauss(st, n_comp, K, sim_state.sigma_s_sq / n_comp);
    x_j = 0;
    if sigma_j_sq > 0
        if sim_state.use_waveform
            j = sim_jammer_waveform(st, n_comp, K, sigma_j_sq, f_j_true, sim_state.fs_hz);
        else
            j = cgauss(st, n_comp, K, sigma_j_sq / n_comp);
        end
        x_j = e_j * j;
        x   = x + x_j;
    end
    x = x + cgauss(st, n_el, K, sim_state.sigma_n_sq);

    % adaptation_tap 'post' = single-path receiver: the monitoring path sits
    % DOWNSTREAM of the notch, so the notch also erases the jammer from the
    % data the DoA / covariance / frequency trackers depend on. This is the
    % blinding failure mode; 'pre' (the recommended architecture) taps upstream
    % and leaves the snapshots untouched. Signal and noise are scaled uniformly
    % rather than spectrally shaped — a power-level approximation consistent
    % with the SINR maths above, and negligible at loss_frac of order 1%.
    if ~isempty(sim_state.notch) && strcmp(sim_state.notch.adaptation_tap, 'post') ...
            && (loss_frac > 0 || h2_jam < 1)
        x = sqrt(1 - loss_frac) * (x - x_j) + sqrt(h2_jam) * x_j;
    end

    obs.snapshots = x;
    if sim_state.use_waveform
        obs.fs_hz = sim_state.fs_hz;
    end
end
end


function [h2_jam, loss_frac, f_notch_used] = notch_state(sim_state, f_notch_cmd, f_j_true)
% Resolve the commanded notch into (jammer power gain, wideband insertion loss).
h2_jam       = 1;
loss_frac    = 0;
f_notch_used = NaN;
if isempty(sim_state.notch) || strcmp(sim_state.notch.mode, 'off')
    return
end
switch sim_state.notch.mode
    case 'oracle'
        f_notch_used = f_j_true;          % perfect knowledge: the upper bound
    case 'adaptive'
        if isempty(f_notch_cmd) || ~isfinite(f_notch_cmd)
            return                        % nothing commanded yet -> not engaged
        end
        f_notch_used = f_notch_cmd;
end
[h2_jam, loss_frac] = sim_notch_response(sim_state.notch, f_notch_used, ...
    f_j_true, sim_state.fs_hz);
end


function x = cgauss(stream, m, n, sigma_sq)
% (m x n) i.i.d. circular complex Gaussian samples with variance sigma_sq.
x = sqrt(sigma_sq / 2) * (randn(stream, m, n) + 1i * randn(stream, m, n));
end
