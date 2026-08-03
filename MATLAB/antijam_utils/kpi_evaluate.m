function kpi = kpi_evaluate(run_log, scenario, antijam_config)
% KPI_EVALUATE  The 5 milestone KPIs from one closed-loop run log.
%
%   kpi = KPI_EVALUATE(run_log, scenario, antijam_config)
%
%   KPI code IS allowed ground truth (scenario.theta_j_deg/phi_j_deg) — the
%   Mode C/S restriction applies to algorithms only.
%
%   Inputs:
%       run_log : struct logged by run_antijam, fields
%           sinr_db        : (1 x T) achieved SINR per step.
%           oracle_sinr_db : (1 x T) perfect-knowledge LCMV SINR per step.
%           W              : (N_el x T) applied weights per step.
%           arm_index      : (1 x T) selected arm (bandit runs; NaN otherwise).
%           grid           : struct with E1 (N_el x N_theta*N_phi), E2 ([] or
%                            same), theta_deg (1 x N_theta), phi_deg
%                            (1 x N_phi), e_s (N_el x n_c) — the full
%                            far-field grid, for pattern-derived KPIs. [P7]
%       scenario       : struct from sim_scenario (truth + events).
%       antijam_config : struct. Required: sinr_min_db.
%
%   Outputs:
%       kpi : struct with fields
%           availability          : fraction of steps with sinr >= sinr_min_db.
%           recovery_steps        : (1 x n_events) steps from each scenario
%                                   event until SINR re-crosses the threshold
%                                   (NaN if never).
%           null_pointing_err_deg : (1 x T) spherical separation from the
%                                   nearest pattern null to (theta_j, phi_j)(t);
%                                   NaN while the jammer is off. [P7: 2-D]
%           peak_gain_penalty_db  : (1 x T) noise-normalized gain toward
%                                   theta_s vs the quiescent MVDR reference
%                                   (jammer-free optimum of the same family).
%           oracle_gap_db         : (1 x T) oracle_sinr_db - sinr_db.
%           oracle_gap_mean_db    : scalar mean of oracle_gap_db.
%           freq_err_hz           : [P10] (1 x T) f_hat_hz - true carrier while
%                                   the jammer is ON; NaN otherwise. Present
%                                   only when the run logged f_hat_hz.
%           freq_rmse_hz          : [P10] scalar RMSE over those steps.
%           notch_gain_db         : [P10] (1 x T) SINR benefit delivered by the
%                                   notch, sinr_db - sinr_no_notch_db. Present
%                                   only when the run logged sinr_no_notch_db.
%                                   Slightly NEGATIVE while the jammer is off —
%                                   that is the insertion loss of a notch held
%                                   through an OFF gap, not an error.
%           notch_gain_on_mean_db : [P10] scalar mean over jammer-ON steps.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6, P7, P10].

thr  = req_field(antijam_config, 'sinr_min_db');
sinr = run_log.sinr_db;
T    = numel(sinr);
grid = run_log.grid;
n_c  = size(grid.e_s, 2);
n_theta = numel(grid.theta_deg);
n_phi   = numel(grid.phi_deg);

% ── 1. Availability ───────────────────────────────────────────────
kpi = struct();
kpi.availability = mean(sinr >= thr);

% ── 2. Recovery per scenario event ────────────────────────────────
n_ev = numel(scenario.events);
kpi.recovery_steps = NaN(1, n_ev);
for i = 1:n_ev
    k0 = find(scenario.t_s >= scenario.events{i}.t_s, 1);
    r  = find(sinr(k0:end) >= thr, 1) - 1;
    if ~isempty(r)
        kpi.recovery_steps(i) = r;
    end
end

% ── 3. Null-pointing error (jammer-on steps only) — [P7] native 2-D ──
kpi.null_pointing_err_deg = NaN(1, T);
for k = 1:T
    if ~scenario.jammer_on(k)
        continue;
    end
    p = reshape(pattern_power(run_log.W(:, k), grid), n_theta, n_phi);
    [mins_theta, mins_phi] = local_minima_2d(p, grid.theta_deg, grid.phi_deg);
    if isempty(mins_theta)
        continue;
    end
    d = angular_separation_deg(mins_theta(:), mins_phi(:), ...
        scenario.theta_j_deg(k), scenario.phi_j_deg(k));
    kpi.null_pointing_err_deg(k) = min(d(:));
end

% ── 4. Peak-gain penalty vs quiescent MVDR reference ──────────────
w_ref = adapt_lcmv(eye(size(grid.E1, 1)), grid.e_s, 0);
g_ref = noise_norm_gain(w_ref, grid.e_s, n_c);
kpi.peak_gain_penalty_db = zeros(1, T);
for k = 1:T
    g = noise_norm_gain(run_log.W(:, k), grid.e_s, n_c);
    kpi.peak_gain_penalty_db(k) = 10 * log10(g / g_ref);
end

% ── 5. Oracle gap ─────────────────────────────────────────────────
kpi.oracle_gap_db      = run_log.oracle_sinr_db - sinr;
kpi.oracle_gap_mean_db = mean(kpi.oracle_gap_db);

% ── 6. [P10] Carrier-estimation error (jammer-on steps only) ──────
% Scored only while the jammer transmits: through an OFF gap the tracker is
% deliberately HOLDING its last estimate (see adapt_freq_update), so scoring it
% against a carrier nobody is radiating would measure the wrong thing.
if isfield(run_log, 'f_hat_hz') && isfield(scenario, 'f_j_hz')
    on = scenario.jammer_on;
    kpi.freq_err_hz     = NaN(1, T);
    kpi.freq_err_hz(on) = run_log.f_hat_hz(on) - scenario.f_j_hz(on);
    scored = ~isnan(kpi.freq_err_hz);
    if any(scored)
        kpi.freq_rmse_hz = sqrt(mean(kpi.freq_err_hz(scored).^2));
    else
        kpi.freq_rmse_hz = NaN;
    end
end

% ── 7. [P10] Notch benefit ────────────────────────────────────────
if isfield(run_log, 'sinr_no_notch_db')
    kpi.notch_gain_db = sinr - run_log.sinr_no_notch_db;
    on = scenario.jammer_on;
    if any(on)
        kpi.notch_gain_on_mean_db = mean(kpi.notch_gain_db(on));
    else
        kpi.notch_gain_on_mean_db = NaN;
    end
end
end


% ────────────────────────── HELPERS ───────────────────────────────

function v = req_field(s, key)
if ~isfield(s, key) || isempty(s.(key))
    error('kpi_evaluate:MissingKey', ...
        'Missing required antijam config key: ''%s''.', key);
end
v = s.(key);
end


function p = pattern_power(w, grid)
% (1 x N_theta*N_phi) grid power pattern, summed over polarization components.
p = abs(w' * grid.E1).^2;
if ~isempty(grid.E2)
    p = p + abs(w' * grid.E2).^2;
end
end


function [mins_theta, mins_phi] = local_minima_2d(p, theta_deg, phi_deg)
% Interior local minima of a (N_theta x N_phi) power grid, 4-neighbor test,
% phi wrapped circularly (no toolbox dependency — base MATLAB only).
n_theta = size(p, 1);
if n_theta < 3
    mins_theta = []; mins_phi = [];
    return
end
center = p(2:end-1, :);
up     = p(1:end-2, :);
down   = p(3:end, :);
if size(p, 2) < 2
    % Degenerate single-phi-column grid (e.g. a 1-D toy pattern): no
    % azimuth neighbor to test against, so theta-only local minima.
    is_min = center < up & center < down;
else
    p_wrap = [p(:, end), p, p(:, 1)];           % circular phi padding
    left   = p_wrap(2:end-1, 1:end-2);
    right  = p_wrap(2:end-1, 3:end);
    is_min = center < up & center < down & center < left & center < right;
end
[ri, ci] = find(is_min);
mins_theta = theta_deg(ri + 1);
mins_phi   = phi_deg(ci);
end


function g = noise_norm_gain(w, e_s, n_c)
% Noise-normalized gain toward theta_s (the SINR-relevant beam quality).
g = (sum(abs(w' * e_s).^2) / n_c) / real(w' * w);
end
