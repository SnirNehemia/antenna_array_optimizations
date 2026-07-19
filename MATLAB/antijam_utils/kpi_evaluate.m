function kpi = kpi_evaluate(run_log, scenario, antijam_config)
% KPI_EVALUATE  The 5 milestone KPIs from one closed-loop run log.
%
%   kpi = KPI_EVALUATE(run_log, scenario, antijam_config)
%
%   KPI code IS allowed ground truth (scenario.theta_j_deg) — the Mode C/S
%   restriction applies to algorithms only.
%
%   Inputs:
%       run_log : struct logged by run_antijam, fields
%           sinr_db        : (1 x T) achieved SINR per step.
%           oracle_sinr_db : (1 x T) perfect-knowledge LCMV SINR per step.
%           W              : (N_el x T) applied weights per step.
%           arm_index      : (1 x T) selected arm (bandit runs; NaN otherwise).
%           cut            : struct with E1 (N_el x N_ang), E2 ([] or same),
%                            angle_deg (1 x N_ang), e_s (N_el x n_c) — the
%                            engine cut, for pattern-derived KPIs.
%       scenario       : struct from sim_scenario (truth + events).
%       antijam_config : struct. Required: sinr_min_db, theta_s_deg.
%
%   Outputs:
%       kpi : struct with fields
%           availability          : fraction of steps with sinr >= sinr_min_db.
%           recovery_steps        : (1 x n_events) steps from each scenario
%                                   event until SINR re-crosses the threshold
%                                   (NaN if never).
%           null_pointing_err_deg : (1 x T) |nearest pattern null - theta_j(t)|;
%                                   NaN while the jammer is off.
%           peak_gain_penalty_db  : (1 x T) noise-normalized gain toward
%                                   theta_s vs the quiescent MVDR reference
%                                   (jammer-free optimum of the same family).
%           oracle_gap_db         : (1 x T) oracle_sinr_db - sinr_db.
%           oracle_gap_mean_db    : scalar mean of oracle_gap_db.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

thr  = req_field(antijam_config, 'sinr_min_db');
sinr = run_log.sinr_db;
T    = numel(sinr);
cut  = run_log.cut;
n_c  = size(cut.e_s, 2);

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

% ── 3. Null-pointing error (jammer-on steps only) ─────────────────
kpi.null_pointing_err_deg = NaN(1, T);
for k = 1:T
    if ~scenario.jammer_on(k)
        continue;
    end
    p = pattern_power(run_log.W(:, k), cut);
    % Local minima of the cut pattern; nearest one to the true jammer angle.
    is_min = [false, p(2:end-1) < p(1:end-2) & p(2:end-1) < p(3:end), false];
    mins = cut.angle_deg(is_min);
    if isempty(mins)
        continue;
    end
    d = abs(mod(mins - scenario.theta_j_deg(k) + 180, 360) - 180);
    kpi.null_pointing_err_deg(k) = min(d);
end

% ── 4. Peak-gain penalty vs quiescent MVDR reference ──────────────
w_ref = adapt_lcmv(eye(size(cut.E1, 1)), cut.e_s, 0);
g_ref = noise_norm_gain(w_ref, cut.e_s, n_c);
kpi.peak_gain_penalty_db = zeros(1, T);
for k = 1:T
    g = noise_norm_gain(run_log.W(:, k), cut.e_s, n_c);
    kpi.peak_gain_penalty_db(k) = 10 * log10(g / g_ref);
end

% ── 5. Oracle gap ─────────────────────────────────────────────────
kpi.oracle_gap_db      = run_log.oracle_sinr_db - sinr;
kpi.oracle_gap_mean_db = mean(kpi.oracle_gap_db);
end


% ────────────────────────── HELPERS ───────────────────────────────

function v = req_field(s, key)
if ~isfield(s, key) || isempty(s.(key))
    error('kpi_evaluate:MissingKey', ...
        'Missing required antijam config key: ''%s''.', key);
end
v = s.(key);
end


function p = pattern_power(w, cut)
% (1 x N_ang) cut power pattern, summed over polarization components.
p = abs(w' * cut.E1).^2;
if ~isempty(cut.E2)
    p = p + abs(w' * cut.E2).^2;
end
end


function g = noise_norm_gain(w, e_s, n_c)
% Noise-normalized gain toward theta_s (the SINR-relevant beam quality).
g = (sum(abs(w' * e_s).^2) / n_c) / real(w' * w);
end
