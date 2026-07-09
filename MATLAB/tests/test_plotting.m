function tests = test_plotting
% TEST_PLOTTING  Smoke tests: each plotting function runs headless and writes
% its expected output file(s). (Plots are visualization, not unit-mathematics.)
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'matlab_utils'));
root = fileparts(fileparts(here));
testCase.TestData.dataDir = fullfile(root, 'data', 'spacing0.6');
outdir = fullfile(tempdir, ['matlab_plot_test_' char(java.util.UUID.randomUUID)]);
mkdir(outdir);
testCase.TestData.outdir = outdir;
end


function teardownOnce(testCase)
if isfolder(testCase.TestData.outdir)
    rmdir(testCase.TestData.outdir, 's');
end
end


function [ep, theta, phi, directives] = load_real(testCase)
patterns = load_element_patterns(testCase.TestData.dataDir);
theta = patterns(1).theta_deg;
phi   = patterns(1).phi_deg;
n     = numel(patterns);
[nt, np_] = size(patterns(1).E_abs);
ep = zeros(n, nt, np_);
for i = 1:n
    ep(i, :, :) = get_component(patterns(i), 'Copol').complex;
end
directives = { ...
    struct('type','peak','theta',30,'phi',0,'theta_width',10,'phi_width',10,'weight',1.0), ...
    struct('type','null','theta',-30,'phi',0,'theta_width',10,'phi_width',10,'weight',10.0) };
end


function oc = plot_config()
oc = struct('plot_cut_type','auto','plot_phi_deg',0,'plot_theta_deg',90, ...
    'plot_dynamic_range_db',40,'plot_equal_area',false, ...
    'save_polar_plot',true,'save_cartesian_plot',true,'save_2d_projection_plot',true, ...
    'save_weight_plots',true,'save_cost_history_plot',true);
end


function test_save_all_plots(testCase)
[ep, theta, phi, directives] = load_real(testCase);
w  = ones(size(ep, 1), 1);
af = compute_array_factor(w, ep);
grid_db = 10 * log10(max(abs(af) .^ 2, 1e-30));

all_cost_histories = {[5; 4; 3], [6; 5; 4; 3.5]};
all_run_labels     = {'Uniform init', 'Random init 1'};

save_all_plots(grid_db, theta, phi, w, directives, ...
    all_cost_histories, 1, all_run_labels, plot_config(), testCase.TestData.outdir);

for f = {'pattern_polar.png','pattern_cartesian.png','pattern_2d.png','weights.png','cost_history.png'}
    verifyTrue(testCase, isfile(fullfile(testCase.TestData.outdir, f{1})), ...
        sprintf('missing %s', f{1}));
end
end


function test_save_pattern_gif(testCase)
[ep, theta, phi, directives] = load_real(testCase);
weights_history = {ones(size(ep,1),1), 0.9*ones(size(ep,1),1), 0.8*ones(size(ep,1),1)};
cost_history    = [-1; -2; -3];
oc = struct('gif_max_frames', 100, 'plot_dynamic_range_db', 40);

save_pattern_gif(weights_history, ep, theta, phi, directives, cost_history, oc, ...
    testCase.TestData.outdir);

verifyTrue(testCase, isfile(fullfile(testCase.TestData.outdir, 'pattern_evolution.gif')));
end


function test_plot_comparison(testCase)
theta = (0:10:180)';
phi   = (0:30:330)';
n_side = 3; d = 0.5;
ep = build_ura_element_patterns(n_side, d, theta, phi);

scenario = struct('label','Smoke scenario','steer_theta_deg',0,'steer_phi_deg',0, ...
    'null_theta_deg',40, ...
    'directives', {{ ...
        struct('type','peak','theta',0,'phi',0,'theta_width',20,'phi_width',360,'weight',1.0), ...
        struct('type','null','theta',40,'phi',0,'theta_width',10,'phi_width',360,'weight',2.0) }});
cfg = struct('n_restarts',1,'max_iterations',50,'cost_tolerance',1e-9, ...
    'd_over_lambda',d,'element_source','synthetic','plot_dynamic_range_db',60);

results = run_scenario(scenario, ep, theta, phi, n_side, d, cfg);
all_scenario_results = { struct('scenario', scenario, 'results', results) };

outpath = fullfile(testCase.TestData.outdir, 'compare.png');
plot_comparison(all_scenario_results, theta, n_side, cfg, outpath);

verifyTrue(testCase, isfile(outpath));
end
