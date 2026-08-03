% RUN_FREQ_NOTCH_DEMO_SCRIPT  Jammer carrier estimation + RF notch [P10] demo.
%
%   Answers the customer's question directly: WHICH FREQUENCY is the unknown
%   jammer on, inside a ~1% band around the array's design frequency, and what
%   does notching it buy on top of the spatial null?
%
%   The demo turns on the [P10] waveform layer (sim.fs_hz), which makes the
%   jammer a CW tone instead of temporally white noise, and runs the real CST
%   array through three cases chosen to bracket the answer:
%
%     D1  OFF-BEAM JAMMER, spatial null works.
%         Establishes that the carrier is estimated accurately, and — the
%         phase's central honest result — that the notch adds almost NOTHING
%         here. The LCMV null has already driven the jammer tens of dB below
%         the noise floor and a notch can only remove what the jammer still
%         contributes. The two mechanisms are complementary coverage, not
%         additive gain.
%
%     D2  MAIN-BEAM JAMMER (theta_j, phi_j = the desired direction).
%         The case the spatial nuller declares OUT OF SCOPE — antijam.guard_deg
%         exists precisely to exclude it, because nulling the main beam would
%         destroy the desired signal. The notch is orthogonal to angle and
%         rescues it outright. This is where the feature earns its place.
%
%     D3  HOPPING CARRIER + ON/OFF jammer.
%         Exercises the tracker's re-lock after a frequency hop and its
%         hold-through-OFF behaviour, and shows the notch covering the
%         covariance tracker's reacquisition lag at each turn-on.
%
%   Each case is run twice — notch "off" (spatial null only) and notch
%   "adaptive" (null + notch commanded by the tracker) — so every number is a
%   like-for-like comparison. D1/D3 additionally report the "oracle" notch
%   (centred on the true carrier) as the upper bound.
%
%   Outputs, in results/freq_notch_demo/<timestamp>/:
%     * freq_waterfall_<case>.png  periodogram vs time with the true carrier,
%                                  the tracked estimate and the notch band;
%                                  the estimation error against the notch
%                                  half-width; and SINR with vs without the
%                                  notch (plot_freq_waterfall);
%     * freq_notch_stats.txt       the printed statistics table.
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P10].

clear; clc;

% ── 0. Paths: reuse the Milestone-1 port + the anti-jam library ─────
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(fileparts(script_dir));           % <repo> = MATLAB/../
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

% ── 1. Config + element patterns ────────────────────────────────────
config = read_config_yaml(fullfile(repo_root, 'config.yaml'));
aj     = config.antijam;
sim_cfg = config.sim;

% Turn ON the [P10] waveform layer for this demo. config.yaml ships with
% sim.fs_hz commented out so the pre-P10 campaign is untouched; 24 MHz is ~1%
% of the 2.4 GHz centre — the customer's stated search window.
sim_cfg.fs_hz      = 24.0e6;
sim_cfg.duration_s = 60.0;

patterns  = load_element_patterns(fullfile(repo_root, config.element_patterns_dir));
theta_deg = patterns(1).theta_deg;
phi_deg   = patterns(1).phi_deg;
[stack1, stack2, pol] = select_polarization_stacks(patterns, config);
fprintf('Array: %d elements, grid %dx%d, polarization %s\n', ...
    size(stack1, 1), numel(theta_deg), numel(phi_deg), pol);
fprintf('Band: f_c = %.3f GHz, captured %.1f MHz (%.2f%% fractional), K = %d\n', ...
    aj.f_center_hz / 1e9, sim_cfg.fs_hz / 1e6, ...
    100 * sim_cfg.fs_hz / aj.f_center_hz, sim_cfg.snapshots_per_step);

% ── 2. Output folder ────────────────────────────────────────────────
timestamp  = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
output_dir = fullfile(repo_root, 'results', 'freq_notch_demo', timestamp);
if ~isfolder(output_dir), mkdir(output_dir); end
fprintf('Output: %s\n\n', output_dir);

% ── 3. The three cases ──────────────────────────────────────────────
% Carrier offsets are deliberately NOT on FFT bin centres (fs/K = 187.5 kHz),
% so the sub-bin interpolation is exercised rather than flattered.
cases = {};
cases{end + 1} = struct( ...
    'name', 'D1_offbeam', ...
    'desc', 'Off-beam jammer: spatial null works, notch adds ~nothing', ...
    'aj_override', struct(), ...
    'scn', struct('id', 'D1', 'motion', 'static', 'power', 'constant', ...
                  'theta_j_deg', 90.0, 'phi_j_deg', 200.0, 'jn_ratio_db', 20.0, ...
                  'freq', 'constant', 'f_j_offset_hz', 3.1e6));
cases{end + 1} = struct( ...
    'name', 'D2_mainbeam', ...
    'desc', 'Main-beam jammer: spatial nulling out of scope, notch rescues it', ...
    ...  % guard_deg = 0 so the scenario generator will place the jammer on the
    ...  % desired direction — normally forbidden, which is the whole point.
    'aj_override', struct('guard_deg', 0.0), ...
    'scn', struct('id', 'D2', 'motion', 'static', 'power', 'constant', ...
                  'theta_j_deg', aj.theta_s_deg, 'phi_j_deg', aj.phi_s_deg, ...
                  'jn_ratio_db', 20.0, ...
                  'freq', 'constant', 'f_j_offset_hz', -5.37e6));
cases{end + 1} = struct( ...
    'name', 'D3_hop_onoff', ...
    'desc', 'Hopping carrier + on/off: re-lock, hold-through-OFF, transient cover', ...
    'aj_override', struct(), ...
    'scn', struct('id', 'D3', 'motion', 'static', 'power', 'onoff', ...
                  'duty_cycle', 0.5, 'toggle_period_s', 10.0, ...
                  'theta_j_deg', 90.0, 'phi_j_deg', 200.0, 'jn_ratio_db', 20.0, ...
                  'freq', 'hop', 'f_j_hop_offsets_hz', [3.1e6, -7.45e6], ...
                  'f_j_hop_period_s', 15.0));

modes = {'off', 'adaptive', 'oracle'};
rows  = {};

for ci = 1:numel(cases)
    c   = cases{ci};
    aj_c = aj;
    f = fieldnames(c.aj_override);
    for i = 1:numel(f), aj_c.(f{i}) = c.aj_override.(f{i}); end
    scn = sim_scenario(c.scn, aj_c, sim_cfg);

    fprintf('=== %s — %s ===\n', c.name, c.desc);
    logs = struct();
    for mi = 1:numel(modes)
        cfg = config;
        cfg.notch = config.notch;
        cfg.notch.mode = modes{mi};
        log = closed_loop_run('lcmv', stack1, stack2, theta_deg, phi_deg, ...
            scn, aj_c, sim_cfg, cfg, []);
        log.oracle_sinr_db = zeros(1, numel(scn.t_s));   % unused by the P10 KPIs
        kpi = kpi_evaluate(log, scn, aj_c);
        logs.(modes{mi}) = log;

        T  = numel(scn.t_s);
        ss = floor(T / 2):T;                              % steady state
        on = scn.jammer_on;
        rows{end + 1} = struct( ...
            'case', c.name, 'notch', modes{mi}, ...
            'sinr_mean_db',  mean(log.sinr_db(ss)), ...
            'availability',  100 * mean(log.sinr_db >= aj_c.sinr_min_db), ...
            'min_sinr_db',   min(log.sinr_db(on)), ...
            'freq_rmse_hz',  kpi_field(kpi, 'freq_rmse_hz'), ...
            'notch_gain_db', kpi_field(kpi, 'notch_gain_on_mean_db')); %#ok<SAGROW>
        fprintf('  notch %-9s SINR %+7.2f dB | avail %5.1f%% | min %+7.2f dB | f RMSE %8s | notch gain %+6.2f dB\n', ...
            modes{mi}, rows{end}.sinr_mean_db, rows{end}.availability, ...
            rows{end}.min_sinr_db, hz_str(rows{end}.freq_rmse_hz), rows{end}.notch_gain_db);
    end
    rescue_db = rows{end - 1}.sinr_mean_db - rows{end - 2}.sinr_mean_db;
    fprintf('  --> adaptive notch vs spatial null only: %+.2f dB\n\n', rescue_db);

    plot_freq_waterfall(logs.adaptive, scn, config.notch, ...
        fullfile(output_dir, ['freq_waterfall_' c.name '.png']));
end

% ── 4. Stats table ──────────────────────────────────────────────────
fid = fopen(fullfile(output_dir, 'freq_notch_stats.txt'), 'w');
hdr = sprintf('%-14s %-9s %12s %8s %11s %12s %12s\n', 'case', 'notch', ...
    'SINR_ss_dB', 'avail_%', 'min_SINR_dB', 'f_RMSE_Hz', 'notch_gain_dB');
fprintf(fid, '[P10] Jammer carrier estimation + RF notch demo\n');
fprintf(fid, 'array %s, pol %s, f_c %.3f GHz, fs %.1f MHz, K %d\n\n', ...
    config.element_patterns_dir, pol, aj.f_center_hz / 1e9, ...
    sim_cfg.fs_hz / 1e6, sim_cfg.snapshots_per_step);
fprintf(fid, '%s', hdr);
for i = 1:numel(rows)
    r = rows{i};
    fprintf(fid, '%-14s %-9s %12.2f %8.1f %11.2f %12s %12.2f\n', ...
        r.case, r.notch, r.sinr_mean_db, r.availability, r.min_sinr_db, ...
        hz_str(r.freq_rmse_hz), r.notch_gain_db);
end
fclose(fid);
fprintf('Wrote %s\n', fullfile(output_dir, 'freq_notch_stats.txt'));


% ────────────────────────── HELPERS ───────────────────────────────

function v = kpi_field(kpi, name)
% KPI fields that only exist when the waveform layer / notch produced them.
if isfield(kpi, name), v = kpi.(name); else, v = NaN; end
end


function s = hz_str(v)
if isnan(v), s = '-'; else, s = sprintf('%.0f', v); end
end
