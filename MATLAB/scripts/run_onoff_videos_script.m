% ==================================================================
%  run_onoff_videos_script.m -- [O] side-by-side videos across the difficulty range
%
%  Renders one MP4 per selected cell: a row of radiation-pattern heatmaps, one
%  per algorithm on a SHARED colour scale, above a single SINR panel carrying
%  every algorithm's trace, the oracle, the operating threshold, and shading
%  over the intervals when the jammer is transmitting.
%
%  The cells are chosen to span the range the campaign measured rather than to
%  flatter the algorithm: one where every method already succeeds, one where
%  the repair earns its keep, one where the reactive path collapses, one on the
%  mirror-degenerate array, and one on an aperture too small to support MUSIC
%  at all (where the honest result is that nothing can be done and the code
%  should degrade rather than fail).
%
%  Part of: Antenna Array Pattern Optimization Tool -- anti-jam milestone [O].
% ==================================================================

function run_onoff_videos_script()

script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(fileparts(script_dir));
addpath(fullfile(script_dir, '..', 'matlab_utils'));
addpath(fullfile(script_dir, '..', 'antijam_utils'));

config = read_config_yaml(fullfile(repo_root, 'config.yaml'));
ONOFF  = struct('fast_lambda', 0.5, 'max_period_s', 60.0, ...
                'lead_frac', 0.15, 'lead_cap_horizons', 2.0);

% id | array | target theta,phi | sep multiple of derived guard | period s | why
CELLS = {
  struct('id','A_easy',      'array','patchs_with_monopoles','th',90,'ph',260,'sep',1.25,'T',25, ...
         'why','comfortable: every method tracks the oracle')
  struct('id','B_repair_pays','array','spacing0.6',          'th',45,'ph',150,'sep',1.25,'T',10, ...
         'why','the repair earns its keep: reactive path is far off the oracle')
  struct('id','C_hard',      'array','spacing0.6_disturbed3','th',60,'ph',260,'sep',1.25,'T',4, ...
         'why','fast toggling, close jammer, perturbed aperture')
  struct('id','D_mirror',    'array','ManyDipoles',          'th',60,'ph',260,'sep',2.0,'T',10, ...
         'why','mirror-degenerate array: the DoA is ambiguous by construction')
  struct('id','E_infeasible','array','patch_back2back',      'th',60,'ph',260,'sep',1.25,'T',10, ...
         'why','2 elements, dual-pol: MUSIC cannot run, must degrade not fail')
};

vid_cfg = struct('max_frames', 160, 'fps', 12, 'dynamic_range_db', 35, ...
                 'format', 'mp4');

out_dir = fullfile(repo_root, 'results', 'onoff_campaign', ...
    [char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss')) '_videos']);
if ~isfolder(out_dir), mkdir(out_dir); end
fprintf('Videos -> %s\n\n', out_dir);

for ic = 1:numel(CELLS)
    C = CELLS{ic};
    fprintf('[%d/%d] %s -- %s\n', ic, numel(CELLS), C.id, C.why);

    patterns = load_element_patterns(fullfile(repo_root, 'data', C.array, filesep));
    pol = pick_polarization(patterns);
    cfg_pol = config; cfg_pol.polarization = pol;
    [s1, s2] = select_polarization_stacks(patterns, cfg_pol);
    th = patterns(1).theta_deg; ph = patterns(1).phi_deg;

    prof = kpi_array_profile(s1, s2, th, ph, C.th, C.ph, config.adapt.diagonal_loading_db);
    sep  = C.sep * prof.guard_deg;
    % Same placement rule the campaign uses, including the phi-cut fallback and
    % the reachability check -- without it a large separation silently clamps to
    % a pole and the label reports an angle the geometry cannot represent.
    [th_j, ph_j, ok] = place_jammer(th, ph, C.th, C.ph, sep);
    if ~ok
        error('run_onoff_videos:UnreachableSeparation', ...
            ['Cell %s asks for %.0f deg separation (%.2f x guard %.0f) on %s, ' ...
             'which is not representable on this grid. Choose a smaller ' ...
             'multiple.'], C.id, sep, C.sep, prof.guard_deg, C.array);
    end

    aj = config.antijam;
    aj.theta_s_deg = C.th; aj.phi_s_deg = C.ph;
    aj.sigma_s_db  = 10.0; aj.jn_ratio_db = 25.0;
    aj.guard_deg   = prof.guard_deg;

    sim_cfg = config.sim;
    sim_cfg.duration_s = 4 * C.T;
    sim_cfg.seed = config.sim.seed;

    scn_cfg = struct('id', 'ONOFF', 'motion', 'static', 'power', 'onoff', ...
        'duty_cycle', 0.5, 'toggle_period_s', C.T, ...
        'theta_j_deg', th_j, 'phi_j_deg', ph_j, 'jn_ratio_db', aj.jn_ratio_db);
    scn = sim_scenario(scn_cfg, aj, sim_cfg);

    cfg_base = config; cfg_base.polarization = pol;
    cfg_base.adapt.predict.doa_stride = 1 + (mean(diff(th)) < 2);
    cfg_fix  = cfg_base;
    cfg_fix.adapt.predict.onoff = ONOFF;

    o = closed_loop_run('oracle',  s1, s2, th, ph, scn, aj, sim_cfg, cfg_base, []);
    l = closed_loop_run('lcmv',    s1, s2, th, ph, scn, aj, sim_cfg, cfg_base, []);
    b = closed_loop_run('predict', s1, s2, th, ph, scn, aj, sim_cfg, cfg_base, []);
    f = closed_loop_run('predict', s1, s2, th, ph, scn, aj, sim_cfg, cfg_fix,  []);
    l.oracle_sinr_db = o.sinr_db;
    b.oracle_sinr_db = o.sinr_db;
    f.oracle_sinr_db = o.sinr_db;

    sc = @(x) 100 * mean((o.sinr_db - x) <= 3.0);
    fprintf(['      %s (%.0f,%.0f) sep %.0f deg (%.2fx guard %.0f), T=%.0f s | ' ...
             'oracle 100  lcmv %.1f  predict %.1f  predict+repair %.1f\n'], ...
        C.array, C.th, C.ph, sep, C.sep, prof.guard_deg, C.T, ...
        sc(l.sinr_db), sc(b.sinr_db), sc(f.sinr_db));

    label = sprintf('%s  %s (%.0f,%.0f)  sep %.0f deg  T=%.0f s', ...
        C.id, C.array, C.th, C.ph, sep, C.T);
    out = fullfile(out_dir, sprintf('%s_%s.mp4', C.id, C.array));
    % The oracle is drawn by the renderer as the dashed reference, taken from
    % run_logs{1}.oracle_sinr_db -- so it must NOT also be passed as a series,
    % or it appears twice in the legend and twice on the plot.
    save_comparison_video({l, b, f}, ...
        {'lcmv', 'predict', 'predict + repair'}, ...
        scn, s1, s2, th, ph, aj, vid_cfg, label, out);
    fprintf('      wrote %s\n\n', out);
end
fprintf('Done. %d videos in %s\n', numel(CELLS), out_dir);
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
