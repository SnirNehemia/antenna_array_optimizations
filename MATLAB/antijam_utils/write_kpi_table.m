function write_kpi_table(runs, aj, output_dir)
% WRITE_KPI_TABLE  Aggregate KPIs per scenario x algorithm; write CSV + txt.
%
%   WRITE_KPI_TABLE(runs, aj, output_dir)
%
%   Aggregation: availability mean, recovery median (all events), null-error
%   median (jammer-on steps), gain-penalty mean, oracle-gap mean. Also prints
%   the txt table to the command window. Shared by run_antijam and
%   run_jammer_demo.
%
%   Inputs:
%       runs       : cell array of run records (scenario_id, algorithm, kpi).
%       aj         : antijam config (sinr_min_db for the header line).
%       output_dir : destination folder (kpi_table.csv / kpi_table.txt).
%
%   Part of: Antenna Array Pattern Optimization Tool — anti-jam milestone [P6].

key = cellfun(@(r) sprintf('%s|%s', r.scenario_id, r.algorithm), runs, ...
              'UniformOutput', false);
[groups, ~, gidx] = unique(key, 'stable');
lines_csv = {['scenario,algorithm,n_runs,availability_mean,recovery_median_steps,' ...
              'null_err_median_deg,gain_penalty_mean_db,oracle_gap_mean_db']};
lines_txt = {sprintf('%-8s %-8s %6s %10s %10s %10s %10s %10s', ...
    'scn', 'alg', 'runs', 'avail', 'recov', 'null_err', 'gain_pen', 'ora_gap')};
for g = 1:numel(groups)
    members = runs(gidx == g);
    parts = strsplit(groups{g}, '|');
    avail = mean(cellfun(@(r) r.kpi.availability, members));
    recs  = cellfun(@(r) {r.kpi.recovery_steps}, members);
    recs  = [recs{:}];
    rec_med = median(recs(~isnan(recs)));
    if isempty(rec_med) || isnan(rec_med), rec_med = NaN; end
    nerr  = cellfun(@(r) {r.kpi.null_pointing_err_deg}, members);
    nerr  = [nerr{:}];
    nerr_med = median(nerr(~isnan(nerr)));
    gpen  = mean(cellfun(@(r) mean(r.kpi.peak_gain_penalty_db), members));
    ogap  = mean(cellfun(@(r) r.kpi.oracle_gap_mean_db, members));
    lines_csv{end + 1} = sprintf('%s,%s,%d,%.4f,%.1f,%.2f,%.2f,%.2f', ...
        parts{1}, parts{2}, numel(members), avail, rec_med, nerr_med, gpen, ogap); %#ok<AGROW>
    lines_txt{end + 1} = sprintf('%-8s %-8s %6d %10.3f %10.1f %10.2f %10.2f %10.2f', ...
        parts{1}, parts{2}, numel(members), avail, rec_med, nerr_med, gpen, ogap); %#ok<AGROW>
end
fid = fopen(fullfile(output_dir, 'kpi_table.csv'), 'w');
fprintf(fid, '%s\n', lines_csv{:});
fclose(fid);
fid = fopen(fullfile(output_dir, 'kpi_table.txt'), 'w');
fprintf(fid, 'Anti-jam KPI table (threshold %.1f dB)\n\n', aj.sinr_min_db);
fprintf(fid, '%s\n', lines_txt{:});
fclose(fid);
fprintf('%s\n', lines_txt{:});
end
