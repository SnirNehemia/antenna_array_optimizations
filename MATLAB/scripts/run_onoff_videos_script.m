% ==================================================================
%  run_onoff_videos_script.m -- [O2] on/off videos, 8 cycles, repair enabled
%
%  Two families:
%
%  A. METHOD COMPARISON (save_comparison_video). One row of radiation patterns,
%     one per algorithm on a SHARED colour scale, above a common SINR panel with
%     the achievable bound, the operating threshold and the ON windows shaded.
%     Cells span the range the campaign measured rather than flattering the
%     algorithm: one where everything already works, one where the repair earns
%     its keep, one hard case, one on the mirror-degenerate array, and one on an
%     aperture too small for MUSIC (where the honest result is that nothing can
%     be done and the code must degrade rather than fail).
%
%  B. AMPLITUDE REGIMES (save_amplitude_grid_video). A 2x2 of weak/strong
%     desired signal against weak/strong jammer, each quadrant showing its
%     pattern above its own SINR trace. Both the colour scale and the SINR axis
%     range are shared across the four, so the comparison shows the true
%     difference in received SINR instead of four independently autoscaled
%     panels that all look equally healthy.
%
%  [O2] Runs are 8 toggle cycles, matching the campaign: period detection needs
%  min_periods = 3 cycles, so the 4-cycle renders of Phase O showed the
%  algorithm mostly in its unlearned state. The on/off repair is now enabled by
%  default in config.yaml, so the 'predict + repair' panel is what ships.
%
%  Part of: Antenna Array Pattern Optimization Tool -- anti-jam milestone [O2].
% ==================================================================

function run_onoff_videos_script(families, cell_idx, out_dir_override)
% FAMILIES  : optional cellstr subset of {'compare','amplitude'} ({} = both).
%             Useful when only one renderer has changed and the other set is
%             still good.
% CELL_IDX  : optional numeric index into the selected family's cell list, so a
%             single video can be rendered in its own process. This machine runs
%             several MATLAB sessions at once and MUSIC's full-grid scan
%             exhausted memory when four 16-element arrays were rendered inside
%             one process; a fresh process per video keeps the footprint flat.

if nargin < 1 || isempty(families), families = {'compare', 'amplitude'}; end
if nargin < 2, cell_idx = []; end
% OUT_DIR_OVERRIDE lets several single-cell invocations write into ONE folder,
% instead of each stamping its own.
if nargin < 3, out_dir_override = ''; end

script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(fileparts(script_dir));
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

config = read_config_yaml(fullfile(repo_root, 'config.yaml'));
CYCLES = 8;

% The shipped repair (config.yaml now carries it); kept explicit so the videos
% do not silently change meaning if the default moves again.
ONOFF = config.adapt.predict.onoff;

SIGMA_S_DB  = 10.0;
JN_RATIO_DB = 25.0;

% ---- A. method-comparison cells -------------------------------------------
CELLS = {
  struct('id','A_easy',       'array','patchs_with_monopoles','th',90,'ph',260,'sep',1.25,'T',25, ...
         'why','comfortable: every method tracks the achievable bound')
  struct('id','B_repair_pays','array','spacing0.6',           'th',45,'ph',150,'sep',1.25,'T',10, ...
         'why','the repair earns its keep: the reactive path is far off the bound')
  struct('id','C_hard',       'array','spacing0.6_disturbed3','th',60,'ph',260,'sep',1.25,'T',4, ...
         'why','fast toggling, close jammer, perturbed aperture')
  struct('id','D_mirror',     'array','ManyDipoles',          'th',60,'ph',260,'sep',2.0, 'T',10, ...
         'why','mirror-degenerate array: the direction estimate is ambiguous by construction')
  struct('id','E_infeasible', 'array','patch_back2back',      'th',60,'ph',260,'sep',1.25,'T',10, ...
         'why','2 elements, dual-pol: MUSIC cannot run, must degrade not fail')
};

% ---- B. amplitude-regime cells (array x geometry) --------------------------
% Chosen to span aperture size, grid resolution and jammer geometry.
AMP_CELLS = {
  struct('id','G1', 'array','spacing0.6',           'th',45,'ph',150,'sep',1.25,'T',10, ...
         'why','16 elements, close jammer -- the hardest geometry in the suite')
  struct('id','G2', 'array','patchs_with_monopoles','th',90,'ph',260,'sep',2.5, 'T',10, ...
         'why','6 elements, very broad beam, well-separated jammer')
  struct('id','G3', 'array','ManyDipoles',          'th',60,'ph',260,'sep',2.0, 'T',10, ...
         'why','20 elements on a 5 deg grid, mirror-degenerate')
  struct('id','G4', 'array','Monopoles',            'th',45,'ph',150,'sep',1.25,'T',10, ...
         'why','14 elements, close jammer')
};
% (weak/strong signal) x (weak/strong jammer), read row-major into the 2x2.
AMP_GRID = [ 0 10;  0 30; 20 10; 20 30];

% Honour a single-cell request so each video can run in its own process.
% NOTE: build these as explicit (possibly empty) index lists. Scaling a range
% by a 0/1 flag looks tidy but `1:n * 0` is empty while `sel * 0` is [0], which
% then indexes element zero and errors.
if any(strcmp(families, 'compare'))
    sel_cmp = 1:numel(CELLS);
    if ~isempty(cell_idx), sel_cmp = cell_idx; end
else
    sel_cmp = [];
end
if any(strcmp(families, 'amplitude'))
    sel_amp = 1:numel(AMP_CELLS);
    if ~isempty(cell_idx), sel_amp = cell_idx; end
else
    sel_amp = [];
end

vid_cfg  = struct('max_frames', 150, 'fps', 12, 'dynamic_range_db', 35, 'format', 'mp4');
grid_cfg = struct('max_frames', 110, 'fps', 10, 'dynamic_range_db', 35, 'format', 'mp4');

if isempty(out_dir_override)
    out_dir = fullfile(repo_root, 'results', 'onoff_campaign', ...
        [char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss')) '_videos8']);
else
    out_dir = out_dir_override;
end
if ~isfolder(out_dir), mkdir(out_dir); end
fprintf('Videos -> %s\n\n', out_dir);

% ══════════════════════════════════════════════════════════════════
%  A. method comparison
% ══════════════════════════════════════════════════════════════════
for ic = sel_cmp
    C = CELLS{ic};
    fprintf('[A %d/%d] %s -- %s\n', ic, numel(CELLS), C.id, C.why);
    [s1, s2, th, ph, pol] = load_array(repo_root, config, C.array);
    prof = kpi_array_profile(s1, s2, th, ph, C.th, C.ph, config.adapt.diagonal_loading_db);
    sep  = C.sep * prof.guard_deg;
    [th_j, ph_j, ok] = place_jammer(th, ph, C.th, C.ph, sep);
    if ~ok
        error('run_onoff_videos:UnreachableSeparation', ...
            'Cell %s asks for %.0f deg separation on %s, not representable here.', ...
            C.id, sep, C.array);
    end

    [aj, sim_cfg, scn] = build_case(config, C, prof, th_j, ph_j, ...
        SIGMA_S_DB, JN_RATIO_DB, CYCLES);
    [cfg_base, cfg_fix] = build_cfgs(config, pol, th, ONOFF);

    o = closed_loop_run('oracle',  s1, s2, th, ph, scn, aj, sim_cfg, cfg_base, []);
    l = closed_loop_run('lcmv',    s1, s2, th, ph, scn, aj, sim_cfg, cfg_base, []);
    b = closed_loop_run('predict', s1, s2, th, ph, scn, aj, sim_cfg, cfg_base, []);
    f = closed_loop_run('predict', s1, s2, th, ph, scn, aj, sim_cfg, cfg_fix,  []);
    l.oracle_sinr_db = o.sinr_db;
    b.oracle_sinr_db = o.sinr_db;
    f.oracle_sinr_db = o.sinr_db;

    sc = @(x) 100 * mean((o.sinr_db - x) <= 3.0);
    fprintf(['      %s (%.0f,%.0f) sep %.0f deg (%.2fx guard %.0f), T=%.0f s, %d cycles | ' ...
             'lcmv %.1f  predict %.1f  predict+repair %.1f\n'], ...
        C.array, C.th, C.ph, sep, C.sep, prof.guard_deg, C.T, CYCLES, ...
        sc(l.sinr_db), sc(b.sinr_db), sc(f.sinr_db));

    label = sprintf('%s  %s (%.0f,%.0f)  sep %.0f deg  T=%.0f s', ...
        C.id, C.array, C.th, C.ph, sep, C.T);
    out = fullfile(out_dir, sprintf('%s_%s.mp4', C.id, C.array));
    % The achievable bound is drawn by the renderer from run_logs{1}, so the
    % oracle must NOT also be passed as a series.
    save_comparison_video({l, b, f}, {'lcmv', 'predict', 'predict + repair'}, ...
        scn, s1, s2, th, ph, aj, vid_cfg, label, out);
    fprintf('      wrote %s\n\n', out);
end

% ══════════════════════════════════════════════════════════════════
%  B. amplitude regimes
% ══════════════════════════════════════════════════════════════════
for ic = sel_amp
    C = AMP_CELLS{ic};
    fprintf('[B %d/%d] %s  %s -- %s\n', ic, numel(AMP_CELLS), C.id, C.array, C.why);
    [s1, s2, th, ph, pol] = load_array(repo_root, config, C.array);
    prof = kpi_array_profile(s1, s2, th, ph, C.th, C.ph, config.adapt.diagonal_loading_db);
    sep  = C.sep * prof.guard_deg;
    [th_j, ph_j, ok] = place_jammer(th, ph, C.th, C.ph, sep);
    if ~ok
        error('run_onoff_videos:UnreachableSeparation', ...
            'Cell %s asks for %.0f deg separation on %s, not representable here.', ...
            C.id, sep, C.array);
    end
    [~, cfg_fix] = build_cfgs(config, pol, th, ONOFF);
    % Coarser MUSIC grid for this family. It compares SINR across amplitude
    % regimes, not angular precision, and the full-grid scan on a 1 deg array
    % is what exhausted memory when several ran in one process.
    cfg_fix.adapt.predict.doa_stride = 2 * cfg_fix.adapt.predict.doa_stride;

    cases = struct('log', {}, 'oracle', {}, 'scn', {}, 'sigma_s_db', {}, 'jn_ratio_db', {});
    txt = '';
    for iq = 1:4
        [aj, sim_cfg, scn] = build_case(config, C, prof, th_j, ph_j, ...
            AMP_GRID(iq, 1), AMP_GRID(iq, 2), CYCLES);
        o = closed_loop_run('oracle',  s1, s2, th, ph, scn, aj, sim_cfg, cfg_fix, []);
        g = closed_loop_run('predict', s1, s2, th, ph, scn, aj, sim_cfg, cfg_fix, []);
        cases(iq) = struct('log', g, 'oracle', o, 'scn', scn, ...
            'sigma_s_db', AMP_GRID(iq, 1), 'jn_ratio_db', AMP_GRID(iq, 2));
        txt = [txt sprintf('  s%.0f/j%.0f: %.1f dB (of %.1f)', ...
            AMP_GRID(iq, 1), AMP_GRID(iq, 2), mean(g.sinr_db), mean(o.sinr_db))]; %#ok<AGROW>
    end
    fprintf('      mean SINR%s\n', txt);

    aj_lbl = config.antijam;
    aj_lbl.theta_s_deg = C.th; aj_lbl.phi_s_deg = C.ph;
    label = sprintf('%s  %s (%.0f,%.0f)  jammer sep %.0f deg  T=%.0f s', ...
        C.id, C.array, C.th, C.ph, sep, C.T);
    out = fullfile(out_dir, sprintf('%s_amplitudes_%s.mp4', C.id, C.array));
    save_amplitude_grid_video(cases, s1, s2, th, ph, aj_lbl, grid_cfg, label, out);
    fprintf('      wrote %s\n\n', out);
end

fprintf('Done. %d comparison + %d amplitude-grid videos in %s\n', ...
    numel(CELLS), numel(AMP_CELLS), out_dir);
end


% ────────────────────────── HELPERS ───────────────────────────────

function [s1, s2, th, ph, pol] = load_array(repo_root, config, name)
patterns = load_element_patterns(fullfile(repo_root, 'data', name, filesep));
pol = pick_polarization(patterns);
cfg_pol = config; cfg_pol.polarization = pol;
[s1, s2] = select_polarization_stacks(patterns, cfg_pol);
th = patterns(1).theta_deg;
ph = patterns(1).phi_deg;
end


function [aj, sim_cfg, scn] = build_case(config, C, prof, th_j, ph_j, ...
                                         sigma_s_db, jn_ratio_db, cycles)
aj = config.antijam;
aj.theta_s_deg = C.th; aj.phi_s_deg = C.ph;
aj.sigma_s_db  = sigma_s_db;
aj.jn_ratio_db = jn_ratio_db;
aj.guard_deg   = prof.guard_deg;          % derived, not the global constant
sim_cfg = config.sim;
sim_cfg.duration_s = cycles * C.T;
sim_cfg.seed = config.sim.seed;
scn_cfg = struct('id', 'ONOFF', 'motion', 'static', 'power', 'onoff', ...
    'duty_cycle', 0.5, 'toggle_period_s', C.T, ...
    'theta_j_deg', th_j, 'phi_j_deg', ph_j, 'jn_ratio_db', jn_ratio_db);
scn = sim_scenario(scn_cfg, aj, sim_cfg);
end


function [cfg_base, cfg_fix] = build_cfgs(config, pol, theta_deg, onoff)
cfg_base = config;
cfg_base.polarization = pol;
cfg_base.adapt.predict.doa_stride = 1 + (mean(diff(theta_deg)) < 2);
% config.yaml now ships the repair ON, so the "un-repaired" reference has to be
% built by REMOVING it rather than by adding it.
if isfield(cfg_base.adapt.predict, 'onoff')
    cfg_base.adapt.predict = rmfield(cfg_base.adapt.predict, 'onoff');
end
cfg_fix = cfg_base;
cfg_fix.adapt.predict.onoff = onoff;
end


function [th_j, ph_j, ok] = place_jammer(theta_deg, phi_deg, th_s, ph_s, sep_deg)
% Identical rule to run_onoff_campaign_script: theta cut first, phi cut as the
% fallback, and a hard "not representable" answer rather than a silent clamp.
ok = true;
th_j = th_s + sep_deg;
ph_j = ph_s;
if th_j > 180 || th_j < 0, th_j = th_s - sep_deg; end
if th_j > 180 || th_j < 0
    th_j = th_s;
    if sind(th_s) < 1e-6, ok = false; return; end
    dphi = sep_deg / sind(th_s);
    if dphi > 180, ok = false; return; end
    ph_j = mod(ph_s + dphi, 360);
end
[it, ip] = nearest_index_2d(theta_deg, phi_deg, th_j, ph_j);
th_j = theta_deg(it);
ph_j = phi_deg(ip);
end


function pol = pick_polarization(patterns)
comps = fieldnames(patterns(1).components);
c1 = stack_component(patterns, comps{1});
p1 = 10 * log10(max(abs(c1(:)).^2) + eps);
p2 = -Inf;
if numel(comps) > 1
    c2 = stack_component(patterns, comps{2});
    p2 = 10 * log10(max(abs(c2(:)).^2) + eps);
end
if (p1 - p2) > 60, pol = comps{1}; else, pol = 'total'; end
end
