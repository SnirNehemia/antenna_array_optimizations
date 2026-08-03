function j = sim_jammer_waveform(stream, n_comp, K, sigma_j_sq, f_bb_hz, fs_hz)
% SIM_JAMMER_WAVEFORM  Jammer amplitude block for one snapshot capture [P10].
%
%   j = SIM_JAMMER_WAVEFORM(stream, n_comp, K, sigma_j_sq, f_bb_hz, fs_hz)
%
%   Returns the (n_comp x K) complex amplitude block that sim_engine_step
%   multiplies by the jammer steering column(s): a continuous-wave TONE at
%   baseband offset f_bb_hz, sampled at fs_hz,
%
%       j_c(t) = sqrt(sigma_j_sq / n_comp) * exp(1i * (2*pi*f_bb*t + phi_c)),
%       t = (0:K-1) / fs_hz,
%
%   with an independent uniform phase phi_c per polarization component, redrawn
%   every step from the simulation's private stream.
%
%   Two properties make this a drop-in replacement for the pre-[P10] i.i.d.
%   complex-Gaussian jammer amplitudes, which is what keeps every P1-P9 gate
%   valid (see antijam_milestone_plan.md P10):
%     1. Per-component power is EXACTLY sigma_j_sq / n_comp in every block (a
%        tone has zero amplitude fluctuation), versus the Gaussian path's 1/K
%        sampling variance. The covariance trackers see strictly better
%        conditioned data, never worse.
%     2. Independent uniform phi_c gives E[j_c j_c'^H] = 0 for c ~= c', so the
%        block-averaged SPATIAL covariance still equals sim_analytic_covariance.
%        Every algorithm consumes only (X*X')/K, so none of them can tell the
%        difference in expectation.
%
%   KNOWN SIMPLIFICATION: a real polarized CW source is coherent ACROSS
%   polarization components (effective jammer rank 1, not 2). The independent
%   per-component phase above deliberately preserves the pre-existing
%   "unpolarized source, independent complex-Gaussian amplitudes per component"
%   convention of sim_engine_init, so sim_analytic_covariance and
%   adapt_music_doa's n_sig = 2*n_comp stay valid. For single-component
%   polarization (n_comp = 1) the question does not arise. Modelling a
%   coherent-polarization jammer is a deferred decision, NOT a bug.
%
%   Inputs:
%       stream     : RandStream handle (sim_state.stream; advances in place).
%       n_comp     : 1 or 2 polarization components.
%       K          : snapshots per step (block length).
%       sigma_j_sq : total jammer power (linear, split equally over components).
%       f_bb_hz    : jammer carrier offset from the array centre frequency [Hz].
%                    Must satisfy |f_bb_hz| <= fs_hz/2 (checked by the caller's
%                    sim_engine_init validation).
%       fs_hz      : snapshot-block sample rate [Hz]. Decoupled from sim.dt_s:
%                    the K columns are a short capture burst (K/fs_hz seconds)
%                    inside each closed-loop step, NOT a span of dt_s.
%
%   Outputs:
%       j : (n_comp x K) complex amplitude block.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P10].

t   = (0:K - 1) / fs_hz;                       % (1 x K) capture-burst time axis
amp = sqrt(sigma_j_sq / n_comp);
phi = 2 * pi * rand(stream, n_comp, 1);        % independent phase per component
j   = amp * exp(1i * (2 * pi * f_bb_hz * t + phi));   % implicit (n_comp x K) expand
end
