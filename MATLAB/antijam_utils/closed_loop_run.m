function log = closed_loop_run(alg, stack1, stack2, theta_deg, phi_deg, scn, aj, sim_cfg, config, codebook)
% CLOSED_LOOP_RUN  One closed-loop run of an algorithm against a scenario.
%
%   log = CLOSED_LOOP_RUN(alg, stack1, stack2, theta_deg, phi_deg, scn, aj, ...
%                         sim_cfg, config, codebook)
%
%   Modes are fixed per algorithm: oracle, lcmv -> Mode C; spsa, bandit ->
%   Mode S. The oracle applies the perfect-knowledge LCMV weights with a
%   one-step lag (R of step k -> w applied at k+1); SPSA warm-starts from the
%   quiescent MVDR beam. Shared by run_antijam and run_jammer_demo.
%
%   [P7]: operates on the full 2-D (theta, phi) far-field grid — no 1-D cut.
%
%   Inputs:
%       alg            : 'oracle' | 'lcmv' | 'spsa' | 'bandit'.
%       stack1, stack2 : (N_el x N_theta x N_phi) far-field stacks (stack2 may
%                        be [] for single-component operation).
%       theta_deg      : (1 x N_theta) elevation grid [deg].
%       phi_deg        : (1 x N_phi) azimuth grid [deg].
%       scn            : struct from sim_scenario.
%       aj             : antijam config section.
%       sim_cfg        : sim config section (seed set per run).
%       config         : full parsed config (adapt / agent sections).
%       codebook       : struct from agent_codebook_build ([] unless alg = bandit).
%
%   [P10]: when the waveform layer is on (sim.fs_hz) a carrier-frequency
%   tracker (adapt_freq_init/update) runs ALONGSIDE the Mode C beamformer and
%   commands the RF notch centre. Like w, the commanded notch is one step
%   stale — this step's capture cannot filter this step's own signal.
%
%   Outputs:
%       log : struct with sinr_db (1 x T), W (N_el x T), arm_index (1 x T,
%             NaN except bandit), grid (E1/E2/theta_deg/phi_deg/e_s), and for
%             bandit runs null_center_deg. Caller adds oracle_sinr_db.
%             [P10] with the waveform layer on, also notch_f_hz and
%             sinr_no_notch_db; with a frequency tracker running, also
%             f_hat_hz, freq_present and freq_snr_db.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6, P7, P10].

mode = 'S';
if any(strcmp(alg, {'oracle', 'lcmv', 'predict'})), mode = 'C'; end

% [P10] RF notch section is optional; sim_engine_init validates it.
notch_cfg = [];
if isfield(config, 'notch') && ~isempty(config.notch)
    notch_cfg = config.notch;
end
st   = sim_engine_init(stack1, stack2, theta_deg, phi_deg, scn, aj, sim_cfg, mode, notch_cfg);
n_el = size(stack1, 1);
T    = numel(scn.t_s);

% [P10] Carrier tracking needs snapshots (Mode C) AND a spectral line to find.
track_freq = strcmp(mode, 'C') && st.use_waveform && ...
    isfield(config.adapt, 'freq') && ~isempty(config.adapt.freq);
if ~isempty(notch_cfg) && strcmp(notch_cfg.mode, 'adaptive') && ~track_freq
    error('closed_loop_run:NoFreqTracker', ...
        ['notch.mode = ''adaptive'' needs a carrier tracker to command it, but ' ...
         'algorithm ''%s'' runs in Mode %s with use_waveform = %d and adapt.freq ' ...
         'present = %d. Use notch.mode ''oracle'' or ''off'' for Mode S, or set ' ...
         'sim.fs_hz and adapt.freq.'], alg, mode, st.use_waveform, ...
        isfield(config.adapt, 'freq'));
end
if track_freq
    fr = adapt_freq_init(config.adapt);
end
f_notch = NaN;   % nothing commanded until the tracker locks

switch alg
    case 'oracle'
        w = adapt_lcmv(eye(n_el), st.e_s, 0);           % quiescent start
    case 'lcmv'
        trk = adapt_tracking_init(config.adapt, st.e_s, n_el);
        w = trk.w;
    case 'predict'
        % P8 MUSIC + on/off-prediction Mode C nuller.
        pr = adapt_predict_init(config.adapt, aj, st.e_s, stack1, stack2, ...
            theta_deg, phi_deg, n_el);
        w = pr.w;
    case 'spsa'
        % Warm start from the quiescent MVDR beam — the uniform start sits
        % ~10-20 dB deeper on the real array and dominates the probe budget.
        w_q = adapt_lcmv(eye(n_el), st.e_s, 0);
        sp = adapt_spsa_init(config.adapt.spsa, w_q, sim_cfg.seed + 5000);
        w = sp.w;
    case 'bandit'
        bd = agent_bandit_init(config.agent, codebook, sim_cfg.seed + 6000);
        w = bd.w;
    otherwise
        error('closed_loop_run:BadAlgorithm', 'Unknown algorithm ''%s''.', alg);
end

log = struct();
log.sinr_db   = zeros(1, T);
log.W         = complex(zeros(n_el, T));
log.arm_index = NaN(1, T);
if strcmp(alg, 'predict')
    % MUSIC/prediction diagnostics (for KPIs and plot_doa_waterfall).
    log.doa_theta_deg = NaN(1, T);
    log.doa_phi_deg   = NaN(1, T);
    log.present       = false(1, T);
    log.predicted_on  = false(1, T);
    log.period_est    = NaN(1, T);
end
if st.use_waveform
    % [P10] notch diagnostics (truth-side, from sim_state.last — never obs).
    log.notch_f_hz       = NaN(1, T);
    log.sinr_no_notch_db = NaN(1, T);
end
if track_freq
    log.f_hat_hz     = NaN(1, T);
    log.freq_present = false(1, T);
    log.freq_snr_db  = NaN(1, T);
    % Periodogram history for plot_freq_waterfall. Kept on a stride derived
    % from T (<= 300 columns) rather than every step: a full-rate history is
    % nfft x T doubles per run, which a Monte Carlo campaign would multiply by
    % scenarios x algorithms x seeds. 300 columns is already finer than the
    % plot's pixel resolution. No config key — nothing here is tunable.
    pspec_stride = max(1, ceil(T / 300));
    pspec_k      = 1:pspec_stride:T;
    log.pspec_k  = pspec_k;
    log.pspec    = [];      % (nfft x numel(pspec_k)), allocated on first estimate
    log.pspec_f_hz = [];
    i_pspec = 0;
end
for k = 1:T
    log.W(:, k) = w;
    [obs, st] = sim_engine_step(st, w, f_notch);
    log.sinr_db(k) = obs.sinr_db;
    if st.use_waveform
        log.notch_f_hz(k)       = st.last.f_notch_hz;
        log.sinr_no_notch_db(k) = st.last.sinr_no_notch_db;
    end
    if track_freq
        [f_notch, fr] = adapt_freq_update(fr, obs);
        log.f_hat_hz(k)     = fr.last.f_hat_hz;
        log.freq_present(k) = fr.last.present;
        log.freq_snr_db(k)  = fr.last.snr_db;
        if any(pspec_k == k)
            if isempty(log.pspec)
                log.pspec      = zeros(numel(fr.last.pspec), numel(pspec_k));
                log.pspec_f_hz = fr.last.f_axis_hz;
            end
            i_pspec = i_pspec + 1;
            log.pspec(:, i_pspec) = fr.last.pspec(:);
        end
    end
    switch alg
        case 'oracle'
            w = adapt_lcmv(sim_analytic_covariance(st), st.e_s, 0);
        case 'lcmv'
            [w, trk] = adapt_tracking_update(trk, obs);
        case 'predict'
            [w, pr] = adapt_predict_update(pr, obs);
            log.doa_theta_deg(k) = pr.last.theta_j_deg;
            log.doa_phi_deg(k)   = pr.last.phi_j_deg;
            log.present(k)       = pr.last.present;
            log.predicted_on(k)  = pr.last.predicted_on;
            log.period_est(k)    = pr.last.period_est;
        case 'spsa'
            [w, sp] = adapt_spsa_update(sp, obs);
        case 'bandit'
            log.arm_index(k) = bd.arm_index;
            [w, bd] = agent_bandit_update(bd, obs);
    end
end
log.grid = struct('E1', st.E1, 'E2', st.E2, 'theta_deg', theta_deg, ...
    'phi_deg', phi_deg, 'e_s', st.e_s);
if strcmp(alg, 'bandit')
    log.null_center_deg = codebook.null_center_deg;   % for the arm heatmap
end
end
