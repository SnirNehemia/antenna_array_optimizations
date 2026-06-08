classdef ManualWeightsTuner < handle
% ══════════════════════════════════════════════════════════════════
% MANUALWEIGHTSTUNER
% Interactive MATLAB uifigure UI for manual antenna array weight tuning
% with live 2-D radiation pattern heatmap and metric feedback.
%
% MATLAB port of scripts/manual_weights.py (Python / tkinter).
% Reuses all compute functions from MATLAB/ unchanged.
%
% Usage:
%   tuner = ManualWeightsTuner('config.yaml');
%   tuner.run();   % blocks until the window is closed
%
% Part of: Antenna Array Pattern Optimization Tool.
% ══════════════════════════════════════════════════════════════════

    % ── CONSTANTS ─────────────────────────────────────────────────
    properties (Constant, Access = private)
        DYNAMIC_RANGE_DB     = 40
        LOG10_EPSILON        = 1e-30
        AMP_SLIDER_MIN       = 0.0
        AMP_SLIDER_MAX       = 2.0
        PHASE_SLIDER_MIN_DEG = -180.0
        PHASE_SLIDER_MAX_DEG =  180.0
        INITIAL_AMPLITUDE    = 1.0
        INITIAL_PHASE_DEG    = 0.0
        DEFAULT_DBI_MIN      = -40.0
        DEFAULT_DBI_MAX      =  10.0
        POLARIZATION_COPOL   = 'copol'
        POLARIZATION_CROSS   = 'cross'
        POLARIZATION_TOTAL   = 'total'
        DISPLAY_RELATIVE     = 'relative'
        DISPLAY_ABSOLUTE     = 'absolute'
        DIRECTIVE_ROW_HEIGHT = 30    % pixels per directive table row
        WEIGHTS_PANEL_HEIGHT = 260   % visible height of scrollable weights area (px)
        RIGHT_PANEL_WIDTH    = 530   % right-column width (px)
        FIG_WIDTH            = 1440  % initial figure width  (px)
        FIG_HEIGHT           = 900   % initial figure height (px)
    end

    % ── PUBLIC PROPERTIES ─────────────────────────────────────────
    properties (Access = public)
        % Set by on_load_config; consumed by the entry-point loop.
        next_config_path = []
    end

    % ── PRIVATE PROPERTIES ────────────────────────────────────────
    properties (Access = private)
        config_path_str
        config
        theta_deg
        phi_deg
        n_elements
        element_patterns_copol
        element_patterns_cross
        weights_complex
        model_name

        active_polarization
        display_mode
        dbi_grid
        last_display_grid
        directive_data
        directive_widget_rows
        syncing_weight_display

        fig
        source_label
        pol_dropdown
        disp_dropdown
        dbi_min_field
        dbi_max_field

        ax
        heatmap_surf
        cbar
        hover_text_h

        status_label

        amp_fields
        amp_sliders
        phase_fields
        phase_sliders

        directives_inner_grid

        label_cost
        label_peak
        label_peak_angle
        label_hpbw
    end

    % ══════════════════════════════════════════════════════════════
    % PUBLIC METHODS
    % ══════════════════════════════════════════════════════════════
    methods (Access = public)

        function obj = ManualWeightsTuner(config_path)
        % Constructor: load config + element patterns, build UI, initial render.
            matlab_dir = fileparts(mfilename('fullpath'));
            if ~contains(path, matlab_dir, 'IgnoreCase', true)
                addpath(matlab_dir);
            end

            obj.config_path_str = char(config_path);
            obj.config          = read_config_yaml(obj.config_path_str);
            obj.validate_config();

            patterns_dir = obj.config.element_patterns_dir;
            if ~isfolder(patterns_dir)
                repo_root    = fileparts(matlab_dir);
                patterns_dir = fullfile(repo_root, patterns_dir);
            end
            patterns = load_element_patterns(patterns_dir);

            obj.theta_deg  = patterns(1).theta_deg(:);
            obj.phi_deg    = patterns(1).phi_deg(:);
            obj.n_elements = numel(patterns);
            obj.model_name = ManualWeightsTuner.dir_basename(patterns_dir);

            n_th = numel(obj.theta_deg);
            n_ph = numel(obj.phi_deg);
            obj.element_patterns_copol = zeros(obj.n_elements, n_th, n_ph);
            obj.element_patterns_cross = zeros(obj.n_elements, n_th, n_ph);
            for k = 1:obj.n_elements
                obj.element_patterns_copol(k, :, :) = patterns(k).E_complex;
                obj.element_patterns_cross(k, :, :) = patterns(k).cross_complex;
            end

            obj.weights_complex        = ones(obj.n_elements, 1, 'like', 1 + 0i);
            obj.active_polarization    = ManualWeightsTuner.POLARIZATION_COPOL;
            obj.display_mode           = ManualWeightsTuner.DISPLAY_RELATIVE;
            obj.directive_data         = {};
            obj.directive_widget_rows  = {};
            obj.syncing_weight_display = false;
            obj.amp_fields             = cell(1, obj.n_elements);
            obj.amp_sliders            = cell(1, obj.n_elements);
            obj.phase_fields           = cell(1, obj.n_elements);
            obj.phase_sliders          = cell(1, obj.n_elements);
            obj.last_display_grid      = [];
            obj.dbi_grid               = [];

            obj.build_ui();

            directives_cfg = {};
            if isfield(obj.config, 'directives')
                directives_cfg = obj.config.directives;
            end
            for k = 1:numel(directives_cfg)
                obj.add_directive_row(directives_cfg{k});
            end

            obj.recompute_and_redraw();
            obj.set_status(sprintf('%d elements loaded  %s  %d directives active', ...
                obj.n_elements, char(183), numel(obj.directive_data)));
        end

        function run(obj)
        % Block until the figure window is closed.
            if isvalid(obj.fig)
                uiwait(obj.fig);
            end
        end

    end

    % ══════════════════════════════════════════════════════════════
    % PRIVATE: UI BUILDERS
    % ══════════════════════════════════════════════════════════════
    methods (Access = private)

        function build_ui(obj)
        % Assemble top-level figure: toolbar row + content (heatmap | controls).
        %
        % KEY LAYOUT RULE: panels that must fill their grid cell must be direct
        % children of a uigridlayout — never children of a plain uipanel — so the
        % grid layout manager can size them. All top-level panels therefore go
        % directly into the appropriate uigridlayout here.
            [~, cfg_name, cfg_ext] = fileparts(obj.config_path_str);
            obj.fig = uifigure( ...
                'Name',            sprintf('Manual Pattern Tuner %s %s', char(8212), [cfg_name cfg_ext]), ...
                'Position',        [80 60 ManualWeightsTuner.FIG_WIDTH ManualWeightsTuner.FIG_HEIGHT], ...
                'Resize',          'on', ...
                'CloseRequestFcn', @(~,~) obj.on_close());

            % Top-level 2-row grid: [toolbar (fixed 56 px); content (fills)].
            top_grid = uigridlayout(obj.fig, [2, 1]);
            top_grid.RowHeight   = {56, '1x'};
            top_grid.ColumnWidth = {'1x'};
            top_grid.Padding     = [0 0 0 0];
            top_grid.RowSpacing  = 0;

            % ── Toolbar ───────────────────────────────────────────
            toolbar_panel = uipanel(top_grid, 'BorderType', 'line');
            toolbar_panel.Layout.Row    = 1;
            toolbar_panel.Layout.Column = 1;
            obj.build_toolbar(toolbar_panel);

            % ── Content: 2-column grid (heatmap | right controls) ──
            % The heatmap panel (build_pattern_panel) and the right-column grid are
            % BOTH direct children of content_grid so the grid manager fills them.
            content_grid = uigridlayout(top_grid, [1, 2]);
            content_grid.Layout.Row    = 2;
            content_grid.Layout.Column = 1;
            content_grid.ColumnWidth   = {'1x', ManualWeightsTuner.RIGHT_PANEL_WIDTH};
            content_grid.Padding       = [4 4 4 4];
            content_grid.ColumnSpacing = 6;

            obj.build_pattern_panel(content_grid);   % places panel at column 1

            % Right column: 3-row grid directly in content_grid.
            right_grid = uigridlayout(content_grid, [3, 1]);
            right_grid.Layout.Row    = 1;
            right_grid.Layout.Column = 2;
            right_grid.RowHeight     = {ManualWeightsTuner.WEIGHTS_PANEL_HEIGHT + 44, 230, 'fit'};
            right_grid.ColumnWidth   = {'1x'};
            right_grid.Padding       = [0 0 0 0];
            right_grid.RowSpacing    = 4;

            obj.build_weights_panel(right_grid);
            obj.build_directives_panel(right_grid);
            obj.build_metrics_panel(right_grid);
        end

        % ── BUILD TOOLBAR ─────────────────────────────────────────

        function build_toolbar(obj, parent)
        % Build toolbar controls:
        %   info block | Load Config | [sep] | Display combo | dBi min/max |
        %   [sep] | Polarisation combo | [sep] | Uniform | Load CSV
        %
        % Labels intentionally have NO explicit BackgroundColor so MATLAB's
        % theme (including dark mode) renders them correctly.
            tg = uigridlayout(parent, [1, 15]);
            tg.ColumnWidth  = {'1x', 100, 10, 60, 88, 58, 62, 58, 62, 10, 80, 68, 10, 105, 120};
            tg.Padding      = [6 4 6 4];
            tg.ColumnSpacing = 2;

            % ── Info block ────────────────────────────────────────
            info_panel = uipanel(tg, 'BorderType', 'none');
            info_panel.Layout.Column = 1;
            info_grid = uigridlayout(info_panel, [2, 1]);
            info_grid.RowHeight  = {'1x', '1x'};
            info_grid.Padding    = [0 0 0 0];
            info_grid.RowSpacing = 0;

            [~, cfg_name, cfg_ext] = fileparts(obj.config_path_str);
            lbl_hdr = uilabel(info_grid, ...
                'Text',       sprintf('Config: %s  %s  Model: %s', [cfg_name cfg_ext], char(183), obj.model_name), ...
                'FontWeight', 'bold', ...
                'FontSize',   11);
            lbl_hdr.Layout.Row = 1;

            obj.source_label = uilabel(info_grid, ...
                'Text',      sprintf('  %s', fileparts(obj.config_path_str)), ...
                'FontSize',  9, ...
                'FontColor', [0.55 0.55 0.55]);
            obj.source_label.Layout.Row = 2;

            % ── Load Config button ────────────────────────────────
            btn_cfg = uibutton(tg, 'Text', 'Load Config...', ...
                'FontSize', 11, ...
                'ButtonPushedFcn', @(~,~) obj.on_load_config());
            btn_cfg.Layout.Column = 2;

            % separators (thin panels)
            mk_sep(tg, 3);

            % ── Display mode ──────────────────────────────────────
            lbl_disp = uilabel(tg, 'Text', 'Display:', 'FontSize', 11, ...
                'HorizontalAlignment', 'right');
            lbl_disp.Layout.Column = 4;

            obj.disp_dropdown = uidropdown(tg, ...
                'Items',           {ManualWeightsTuner.DISPLAY_RELATIVE, ManualWeightsTuner.DISPLAY_ABSOLUTE}, ...
                'Value',           ManualWeightsTuner.DISPLAY_RELATIVE, ...
                'FontSize',        11, ...
                'ValueChangedFcn', @(~,~) obj.on_display_mode_change());
            obj.disp_dropdown.Layout.Column = 5;

            lbl_min = uilabel(tg, 'Text', 'min (dBi):', 'FontSize', 11, ...
                'HorizontalAlignment', 'right');
            lbl_min.Layout.Column = 6;

            obj.dbi_min_field = uieditfield(tg, 'numeric', ...
                'Value',           ManualWeightsTuner.DEFAULT_DBI_MIN, ...
                'FontSize',        11, ...
                'ValueChangedFcn', @(~,~) obj.on_display_mode_change());
            obj.dbi_min_field.Layout.Column = 7;

            lbl_max = uilabel(tg, 'Text', 'max (dBi):', 'FontSize', 11, ...
                'HorizontalAlignment', 'right');
            lbl_max.Layout.Column = 8;

            obj.dbi_max_field = uieditfield(tg, 'numeric', ...
                'Value',           ManualWeightsTuner.DEFAULT_DBI_MAX, ...
                'FontSize',        11, ...
                'ValueChangedFcn', @(~,~) obj.on_display_mode_change());
            obj.dbi_max_field.Layout.Column = 9;

            mk_sep(tg, 10);

            % ── Polarisation ──────────────────────────────────────
            lbl_pol = uilabel(tg, 'Text', 'Polarisation:', 'FontSize', 11, ...
                'HorizontalAlignment', 'right');
            lbl_pol.Layout.Column = 11;

            obj.pol_dropdown = uidropdown(tg, ...
                'Items',           {ManualWeightsTuner.POLARIZATION_COPOL, ...
                                    ManualWeightsTuner.POLARIZATION_CROSS, ...
                                    ManualWeightsTuner.POLARIZATION_TOTAL}, ...
                'Value',           ManualWeightsTuner.POLARIZATION_COPOL, ...
                'FontSize',        11, ...
                'ValueChangedFcn', @(~,~) obj.on_polarization_change());
            obj.pol_dropdown.Layout.Column = 12;

            mk_sep(tg, 13);

            btn_uni = uibutton(tg, 'Text', 'Uniform Weights', ...
                'FontSize', 11, ...
                'ButtonPushedFcn', @(~,~) obj.on_uniform());
            btn_uni.Layout.Column = 14;

            btn_csv = uibutton(tg, 'Text', 'Load Weights CSV', ...
                'FontSize', 11, ...
                'ButtonPushedFcn', @(~,~) obj.on_load_csv());
            btn_csv.Layout.Column = 15;
        end

        % ── BUILD PATTERN PANEL ───────────────────────────────────

        function build_pattern_panel(obj, parent_grid)
        % Embed pcolor heatmap + status bar. Panel placed at column 1 of parent_grid
        % so the uigridlayout manager fills it to its cell.
            lf = uipanel(parent_grid, 'Title', '2-D Radiation Pattern');
            lf.Layout.Row    = 1;
            lf.Layout.Column = 1;

            lf_grid = uigridlayout(lf, [2, 1]);
            lf_grid.RowHeight   = {'1x', 24};
            lf_grid.ColumnWidth = {'1x'};
            lf_grid.Padding     = [2 2 2 2];
            lf_grid.RowSpacing  = 2;

            obj.ax = uiaxes(lf_grid);
            obj.ax.Layout.Row    = 1;
            obj.ax.Layout.Column = 1;
            obj.ax.Color         = [0 0 0];
            obj.ax.XLabel.String = 'Azimuth \phi ({\circ})';
            obj.ax.YLabel.String = 'Elevation \theta ({\circ})';
            obj.ax.Title.String  = sprintf('Array Factor %s %s', char(8212), obj.model_name);
            obj.ax.Title.FontSize = 11;
            obj.ax.XLim = [0 360];
            obj.ax.YLim = [0 180];
            obj.ax.YDir = 'reverse';
            grid(obj.ax, 'on');
            obj.ax.GridAlpha     = 0.20;
            obj.ax.GridColor     = [1 1 1];
            obj.ax.GridLineWidth = 0.5;
            obj.ax.FontSize      = 10;

            n_theta  = numel(obj.theta_deg);
            n_phi    = numel(obj.phi_deg);
            init_grid = zeros(n_theta, n_phi) - ManualWeightsTuner.DYNAMIC_RANGE_DB;
            hold(obj.ax, 'on');
            obj.heatmap_surf = pcolor(obj.ax, obj.phi_deg, obj.theta_deg, init_grid);
            obj.heatmap_surf.EdgeColor = 'none';
            colormap(obj.ax, 'jet');
            obj.cbar = colorbar(obj.ax);
            obj.cbar.Label.String = 'dB (normalized to peak)';
            obj.cbar.Label.FontSize = 10;
            obj.cbar.FontSize = 10;
            obj.cbar.Ticks = -ManualWeightsTuner.DYNAMIC_RANGE_DB:10:0;
            clim(obj.ax, [-ManualWeightsTuner.DYNAMIC_RANGE_DB, 0]);

            obj.hover_text_h = text(obj.ax, 0.01, 0.98, '', ...
                'Units',               'normalized', ...
                'VerticalAlignment',   'top', ...
                'HorizontalAlignment', 'left', ...
                'FontSize',            10, ...
                'Color',               [1 1 1], ...
                'BackgroundColor',     [0 0 0], ...
                'Visible',             'off');

            obj.fig.WindowButtonMotionFcn = @(~,~) obj.on_mouse_move();

            obj.status_label = uilabel(lf_grid, ...
                'Text',                '', ...
                'FontSize',            10, ...
                'HorizontalAlignment', 'left');
            obj.status_label.Layout.Row    = 2;
            obj.status_label.Layout.Column = 1;
        end

        % ── BUILD WEIGHTS PANEL ───────────────────────────────────

        function build_weights_panel(obj, parent)
        % Scrollable panel: one row per element (index | amp field+slider | phase field+slider | Solo).
            outer = uipanel(parent, 'Title', 'Element Weights', 'Scrollable', 'on');
            outer.Layout.Row    = 1;
            outer.Layout.Column = 1;

            n_rows = obj.n_elements + 1;   % +1 for header
            w_grid = uigridlayout(outer, [n_rows, 6]);
            w_grid.ColumnWidth   = {40, 62, '1x', 66, '1x', 46};
            w_grid.RowHeight     = repmat({26}, 1, n_rows);
            w_grid.Padding       = [4 4 4 4];
            w_grid.RowSpacing    = 2;
            w_grid.ColumnSpacing = 3;

            % Header row.
            headers = {'Elem', 'Amplitude', '', 'Phase (°)', '', ''};
            for c = 1:6
                lbl = uilabel(w_grid, 'Text', headers{c}, ...
                    'FontWeight', 'bold', 'FontSize', 10, ...
                    'HorizontalAlignment', 'center');
                lbl.Layout.Row    = 1;
                lbl.Layout.Column = c;
            end

            for n = 1:obj.n_elements
                r = n + 1;

                lbl_idx = uilabel(w_grid, 'Text', sprintf('%3d', n - 1), ...
                    'FontSize', 10, 'HorizontalAlignment', 'right');
                lbl_idx.Layout.Row    = r;
                lbl_idx.Layout.Column = 1;

                amp_ef = uieditfield(w_grid, 'numeric', ...
                    'Value',           ManualWeightsTuner.INITIAL_AMPLITUDE, ...
                    'Limits',          [0 Inf], ...
                    'FontSize',        10, ...
                    'ValueChangedFcn', @(src,~) obj.on_amp_entry_changed(n, src.Value));
                amp_ef.Layout.Row    = r;
                amp_ef.Layout.Column = 2;

                amp_sl = uislider(w_grid, ...
                    'Limits',           [ManualWeightsTuner.AMP_SLIDER_MIN, ManualWeightsTuner.AMP_SLIDER_MAX], ...
                    'Value',            ManualWeightsTuner.INITIAL_AMPLITUDE, ...
                    'MajorTicks',       [], 'MinorTicks', [], ...
                    'ValueChangedFcn',  @(src,~) obj.on_amp_slider_changed(n, src.Value), ...
                    'ValueChangingFcn', @(src,evt) obj.on_amp_slider_changing(n, evt.Value));
                amp_sl.Layout.Row    = r;
                amp_sl.Layout.Column = 3;

                phase_ef = uieditfield(w_grid, 'numeric', ...
                    'Value',           ManualWeightsTuner.INITIAL_PHASE_DEG, ...
                    'FontSize',        10, ...
                    'ValueChangedFcn', @(src,~) obj.on_phase_entry_changed(n, src.Value));
                phase_ef.Layout.Row    = r;
                phase_ef.Layout.Column = 4;

                phase_sl = uislider(w_grid, ...
                    'Limits',           [ManualWeightsTuner.PHASE_SLIDER_MIN_DEG, ManualWeightsTuner.PHASE_SLIDER_MAX_DEG], ...
                    'Value',            ManualWeightsTuner.INITIAL_PHASE_DEG, ...
                    'MajorTicks',       [], 'MinorTicks', [], ...
                    'ValueChangedFcn',  @(src,~) obj.on_phase_slider_changed(n, src.Value), ...
                    'ValueChangingFcn', @(src,evt) obj.on_phase_slider_changing(n, evt.Value));
                phase_sl.Layout.Row    = r;
                phase_sl.Layout.Column = 5;

                solo_btn = uibutton(w_grid, 'Text', 'Solo', ...
                    'FontSize',        10, ...
                    'ButtonPushedFcn', @(~,~) obj.on_solo_element(n));
                solo_btn.Layout.Row    = r;
                solo_btn.Layout.Column = 6;

                obj.amp_fields{n}    = amp_ef;
                obj.amp_sliders{n}   = amp_sl;
                obj.phase_fields{n}  = phase_ef;
                obj.phase_sliders{n} = phase_sl;
            end
        end

        % ── BUILD DIRECTIVES PANEL ────────────────────────────────

        function build_directives_panel(obj, parent)
        % Panel with dynamic add/remove rows.
            dir_outer = uipanel(parent, 'Title', 'Directives');
            dir_outer.Layout.Row    = 2;
            dir_outer.Layout.Column = 1;

            outer_grid = uigridlayout(dir_outer, [3, 1]);
            outer_grid.RowHeight   = {24, '1x', 32};
            outer_grid.ColumnWidth = {'1x'};
            outer_grid.Padding     = [2 2 2 2];
            outer_grid.RowSpacing  = 2;

            % Column header row.
            th = char(952); ph = char(966);
            hdr_g = uigridlayout(outer_grid, [1, 8]);
            hdr_g.Layout.Row    = 1;
            hdr_g.Layout.Column = 1;
            hdr_g.ColumnWidth   = {58, 46, 46, 46, 46, 36, 30, '1x'};
            hdr_g.Padding       = [0 0 0 0];
            hdr_g.ColumnSpacing = 2;
            hdr_names = {'Type', [th,'(°)'], [ph,'(°)'], [th,'W'], [ph,'W'], 'Wt', '', 'Result'};
            for c = 1:8
                hl = uilabel(hdr_g, 'Text', hdr_names{c}, ...
                    'FontWeight', 'bold', 'FontSize', 10, ...
                    'HorizontalAlignment', 'center');
                hl.Layout.Column = c;
            end

            % Scrollable container for dynamic rows.
            dir_scroll = uipanel(outer_grid, 'BorderType', 'none', 'Scrollable', 'on');
            dir_scroll.Layout.Row    = 2;
            dir_scroll.Layout.Column = 1;

            obj.directives_inner_grid = uigridlayout(dir_scroll, [1, 1]);
            obj.directives_inner_grid.ColumnWidth = {'1x'};
            obj.directives_inner_grid.RowHeight   = {ManualWeightsTuner.DIRECTIVE_ROW_HEIGHT};
            obj.directives_inner_grid.Padding     = [0 0 0 0];
            obj.directives_inner_grid.RowSpacing  = 2;

            add_btn = uibutton(outer_grid, 'Text', '＋  Add Directive', ...
                'FontSize',        11, ...
                'ButtonPushedFcn', @(~,~) obj.add_directive_row([]));
            add_btn.Layout.Row    = 3;
            add_btn.Layout.Column = 1;
        end

        % ── BUILD METRICS PANEL ───────────────────────────────────

        function build_metrics_panel(obj, parent)
        % Four fixed metric rows.
            met_panel = uipanel(parent, 'Title', 'Metrics');
            met_panel.Layout.Row    = 3;
            met_panel.Layout.Column = 1;

            met_grid = uigridlayout(met_panel, [4, 2]);
            met_grid.ColumnWidth = {100, '1x'};
            met_grid.RowHeight   = repmat({22}, 1, 4);
            met_grid.Padding     = [6 6 6 6];
            met_grid.RowSpacing  = 3;

            metric_titles = {'Total J:', 'Global peak:', 'Peak angle:', '3 dB HPBW:'};
            val_handles   = cell(1, 4);
            for r = 1:4
                tl = uilabel(met_grid, 'Text', metric_titles{r}, ...
                    'FontSize', 11, 'HorizontalAlignment', 'right');
                tl.Layout.Row    = r;
                tl.Layout.Column = 1;
                val_handles{r} = uilabel(met_grid, 'Text', char(8212), ...
                    'FontSize', 11, 'FontName', 'Courier New');
                val_handles{r}.Layout.Row    = r;
                val_handles{r}.Layout.Column = 2;
            end
            obj.label_cost       = val_handles{1};
            obj.label_peak       = val_handles{2};
            obj.label_peak_angle = val_handles{3};
            obj.label_hpbw       = val_handles{4};
        end

    end % UI builder methods

    % ══════════════════════════════════════════════════════════════
    % PRIVATE: DIRECTIVE TABLE MANAGEMENT
    % ══════════════════════════════════════════════════════════════
    methods (Access = private)

        function add_directive_row(obj, directive_dict)
        % Append a directive row, using directive_dict fields or defaults.
            D = struct('type', 'peak', 'theta', 0.0, 'phi', 0.0, ...
                       'theta_width', 5.0, 'phi_width', 5.0, 'weight', 1.0);
            if ~isempty(directive_dict) && isstruct(directive_dict)
                d = directive_dict;
                sym_w = [];
                if isfield(d, 'width'),       sym_w         = d.width;       end
                if isfield(d, 'type'),        D.type        = d.type;        end
                if isfield(d, 'theta'),       D.theta       = d.theta;       end
                if isfield(d, 'phi'),         D.phi         = d.phi;         end
                if isfield(d, 'theta_width'), D.theta_width = d.theta_width;
                elseif ~isempty(sym_w),       D.theta_width = sym_w;         end
                if isfield(d, 'phi_width'),   D.phi_width   = d.phi_width;
                elseif ~isempty(sym_w),       D.phi_width   = sym_w;         end
                if isfield(d, 'weight'),      D.weight      = d.weight;      end
            end

            k = numel(obj.directive_data) + 1;
            if k == 1
                obj.directives_inner_grid.RowHeight = ...
                    {ManualWeightsTuner.DIRECTIVE_ROW_HEIGHT};
            else
                obj.directives_inner_grid.RowHeight{end + 1} = ...
                    ManualWeightsTuner.DIRECTIVE_ROW_HEIGHT;
            end

            rw = obj.make_directive_row_widgets(k, D);
            obj.directive_data{k}        = D;
            obj.directive_widget_rows{k} = rw;

            obj.on_directive_change();
            obj.update_status_count();
        end

        function rw = make_directive_row_widgets(obj, row_k, D)
        % Create all widgets for directive row row_k from data struct D.
            row_grid = uigridlayout(obj.directives_inner_grid, [1, 8]);
            row_grid.Layout.Row    = row_k;
            row_grid.Layout.Column = 1;
            row_grid.ColumnWidth   = {58, 46, 46, 46, 46, 36, 30, '1x'};
            row_grid.Padding       = [0 1 0 1];
            row_grid.ColumnSpacing = 2;

            type_dd = uidropdown(row_grid, ...
                'Items',           {'peak', 'null'}, ...
                'Value',           D.type, ...
                'FontSize',        10, ...
                'ValueChangedFcn', @(src,~) obj.on_directive_field_changed(row_k, 'type', src.Value));
            type_dd.Layout.Column = 1;

            field_names = {'theta', 'phi', 'theta_width', 'phi_width', 'weight'};
            field_vals  = [D.theta, D.phi, D.theta_width, D.phi_width, D.weight];
            efs = cell(1, 5);
            for f = 1:5
                ef = uieditfield(row_grid, 'numeric', ...
                    'Value',           field_vals(f), ...
                    'FontSize',        10, ...
                    'ValueChangedFcn', @(src,~) obj.on_directive_field_changed( ...
                                                    row_k, field_names{f}, src.Value));
                ef.Layout.Column = f + 1;
                efs{f} = ef;
            end

            rm_btn = uibutton(row_grid, 'Text', char(215), ...
                'FontSize',        12, ...
                'ButtonPushedFcn', @(~,~) obj.remove_directive_row(row_k));
            rm_btn.Layout.Column = 7;

            metric_lbl = uilabel(row_grid, 'Text', char(8212), ...
                'FontSize', 10, 'FontName', 'Courier New');
            metric_lbl.Layout.Column = 8;

            rw.row_grid   = row_grid;
            rw.type_dd    = type_dd;
            rw.theta_ef   = efs{1};
            rw.phi_ef     = efs{2};
            rw.tw_ef      = efs{3};
            rw.pw_ef      = efs{4};
            rw.weight_ef  = efs{5};
            rw.metric_lbl = metric_lbl;
        end

        function on_directive_field_changed(obj, row_k, field_name, value)
            if row_k <= numel(obj.directive_data)
                obj.directive_data{row_k}.(field_name) = value;
            end
            obj.on_directive_change();
        end

        function remove_directive_row(obj, row_k)
        % Delete all row widgets, remove data entry, rebuild remaining rows.
            if row_k > numel(obj.directive_data), return; end

            for k = 1:numel(obj.directive_widget_rows)
                if isvalid(obj.directive_widget_rows{k}.row_grid)
                    delete(obj.directive_widget_rows{k}.row_grid);
                end
            end
            obj.directive_widget_rows = {};

            obj.directive_data(row_k) = [];
            n_remain = numel(obj.directive_data);

            if n_remain == 0
                obj.directives_inner_grid.RowHeight = ...
                    {ManualWeightsTuner.DIRECTIVE_ROW_HEIGHT};
            else
                obj.directives_inner_grid.RowHeight = ...
                    repmat({ManualWeightsTuner.DIRECTIVE_ROW_HEIGHT}, 1, n_remain);
                for k = 1:n_remain
                    obj.directive_widget_rows{k} = ...
                        obj.make_directive_row_widgets(k, obj.directive_data{k});
                end
            end

            obj.on_directive_change();
            obj.update_status_count();
        end

    end % directive management methods

    % ══════════════════════════════════════════════════════════════
    % PRIVATE: EVENT HANDLERS
    % ══════════════════════════════════════════════════════════════
    methods (Access = private)

        function on_amp_entry_changed(obj, elem_idx, value)
            if obj.syncing_weight_display, return; end
            amplitude_linear = max(0.0, value);
            phase_deg = angle(obj.weights_complex(elem_idx)) * 180 / pi;
            obj.weights_complex(elem_idx) = amplitude_linear * exp(1i * deg2rad(phase_deg));
            obj.syncing_weight_display = true;
            obj.amp_sliders{elem_idx}.Value = min(amplitude_linear, ManualWeightsTuner.AMP_SLIDER_MAX);
            obj.syncing_weight_display = false;
            obj.recompute_and_redraw();
        end

        function on_amp_slider_changed(obj, elem_idx, value)
            if obj.syncing_weight_display, return; end
            obj.apply_amp_change(elem_idx, value);
        end

        function on_amp_slider_changing(obj, elem_idx, value)
            if obj.syncing_weight_display, return; end
            obj.apply_amp_change(elem_idx, value);
        end

        function apply_amp_change(obj, elem_idx, value)
            amplitude_linear = max(0.0, value);
            phase_deg = angle(obj.weights_complex(elem_idx)) * 180 / pi;
            obj.weights_complex(elem_idx) = amplitude_linear * exp(1i * deg2rad(phase_deg));
            obj.syncing_weight_display = true;
            obj.amp_fields{elem_idx}.Value = amplitude_linear;
            obj.syncing_weight_display = false;
            obj.recompute_and_redraw();
        end

        function on_phase_entry_changed(obj, elem_idx, value)
            if obj.syncing_weight_display, return; end
            phase_deg = value;
            amplitude_linear = abs(obj.weights_complex(elem_idx));
            obj.weights_complex(elem_idx) = amplitude_linear * exp(1i * deg2rad(phase_deg));
            obj.syncing_weight_display = true;
            clamped = max(ManualWeightsTuner.PHASE_SLIDER_MIN_DEG, ...
                      min(ManualWeightsTuner.PHASE_SLIDER_MAX_DEG, phase_deg));
            obj.phase_sliders{elem_idx}.Value = clamped;
            obj.syncing_weight_display = false;
            obj.recompute_and_redraw();
        end

        function on_phase_slider_changed(obj, elem_idx, value)
            if obj.syncing_weight_display, return; end
            obj.apply_phase_change(elem_idx, value);
        end

        function on_phase_slider_changing(obj, elem_idx, value)
            if obj.syncing_weight_display, return; end
            obj.apply_phase_change(elem_idx, value);
        end

        function apply_phase_change(obj, elem_idx, value)
            phase_deg = value;
            amplitude_linear = abs(obj.weights_complex(elem_idx));
            obj.weights_complex(elem_idx) = amplitude_linear * exp(1i * deg2rad(phase_deg));
            obj.syncing_weight_display = true;
            obj.phase_fields{elem_idx}.Value = phase_deg;
            obj.syncing_weight_display = false;
            obj.recompute_and_redraw();
        end

        function on_solo_element(obj, elem_idx)
            obj.weights_complex = zeros(obj.n_elements, 1, 'like', 1 + 0i);
            obj.weights_complex(elem_idx) = 1.0;
            obj.sync_entries_to_weights();
            obj.recompute_and_redraw();
        end

        function on_uniform(obj)
            obj.weights_complex = ones(obj.n_elements, 1, 'like', 1 + 0i);
            obj.sync_entries_to_weights();
            obj.recompute_and_redraw();
        end

        function on_load_csv(obj)
            [fname, fpath] = uigetfile( ...
                {'*.csv', 'CSV files (*.csv)'; '*.*', 'All files (*.*)'}, ...
                'Load Weights CSV');
            if isequal(fname, 0), return; end
            filepath = fullfile(fpath, fname);
            try
                loaded_weights = obj.parse_weights_csv(filepath);
            catch exc
                uialert(obj.fig, exc.message, 'Load Error');
                return;
            end
            obj.weights_complex = loaded_weights;
            obj.sync_entries_to_weights();
            obj.recompute_and_redraw();
            obj.set_status(sprintf('Weights loaded from: %s', filepath));
            config_dir  = fileparts(obj.config_path_str);
            weights_dir = fileparts(filepath);
            if strcmp(config_dir, weights_dir)
                obj.source_label.Text = sprintf('  %s', config_dir);
            else
                obj.source_label.Text = ...
                    sprintf('  Config: %s   Weights: %s', config_dir, weights_dir);
            end
        end

        function on_load_config(obj)
            config_dir = fileparts(obj.config_path_str);
            [fname, fpath] = uigetfile( ...
                {'*.yaml;*.yml', 'YAML files (*.yaml, *.yml)'; '*.*', 'All files (*.*)'}, ...
                'Load Config YAML', config_dir);
            if isequal(fname, 0), return; end
            obj.next_config_path = fullfile(fpath, fname);
            delete(obj.fig);
        end

        function on_directive_change(obj)
            obj.recompute_and_redraw();
        end

        function on_polarization_change(obj)
            obj.active_polarization = obj.pol_dropdown.Value;
            obj.recompute_and_redraw();
        end

        function on_display_mode_change(obj)
            obj.display_mode = obj.disp_dropdown.Value;
            if strcmp(obj.display_mode, ManualWeightsTuner.DISPLAY_ABSOLUTE)
                obj.cbar.Label.String = 'dBi (absolute directivity)';
            else
                obj.cbar.Label.String = 'dB (normalized to peak)';
                obj.cbar.Ticks = -ManualWeightsTuner.DYNAMIC_RANGE_DB:10:0;
            end
            obj.recompute_and_redraw();
        end

        function on_mouse_move(obj)
        % Update cursor annotation with angle and pattern value at cursor.
            if isempty(obj.last_display_grid) || ~isvalid(obj.ax), return; end
            cp           = obj.ax.CurrentPoint;
            phi_cursor   = cp(1, 1);
            theta_cursor = cp(1, 2);
            xl = obj.ax.XLim;  yl = obj.ax.YLim;
            if phi_cursor   < xl(1) || phi_cursor   > xl(2) || ...
               theta_cursor < yl(1) || theta_cursor > yl(2)
                obj.hover_text_h.Visible = 'off';
                return;
            end
            [~, phi_idx]   = min(abs(obj.phi_deg   - phi_cursor));
            [~, theta_idx] = min(abs(obj.theta_deg - theta_cursor));
            val = obj.last_display_grid(theta_idx, phi_idx);
            if strcmp(obj.display_mode, ManualWeightsTuner.DISPLAY_ABSOLUTE)
                unit_str = 'dBi';
            else
                unit_str = 'dB';
            end
            obj.hover_text_h.String = sprintf( ...
                '\\theta = %.1f%s   \\phi = %.1f%s   D = %+.2f %s', ...
                obj.theta_deg(theta_idx), char(176), obj.phi_deg(phi_idx), char(176), val, unit_str);
            obj.hover_text_h.Visible = 'on';
        end

        function on_close(obj)
            if isvalid(obj.fig)
                uiresume(obj.fig);
                delete(obj.fig);
            end
        end

    end % event handler methods

    % ══════════════════════════════════════════════════════════════
    % PRIVATE: RECOMPUTE AND REDRAW
    % ══════════════════════════════════════════════════════════════
    methods (Access = private)

        function recompute_and_redraw(obj)
        % Compute array factor + metrics, refresh all displays.
            w_power = sum(abs(obj.weights_complex) .^ 2);
            w_norm  = obj.weights_complex / sqrt(max(w_power, ManualWeightsTuner.LOG10_EPSILON));

            af_copol = compute_array_factor(w_norm, obj.element_patterns_copol);
            af_cross = compute_array_factor(w_norm, obj.element_patterns_cross);

            % Spherical total radiated power (CST partial-directivity normaliser).
            theta_rad = deg2rad(obj.theta_deg(:));
            if numel(theta_rad) > 1, dtheta_rad = mean(diff(theta_rad)); else, dtheta_rad = pi; end
            if numel(obj.phi_deg) > 1
                dphi_rad = deg2rad(mean(diff(obj.phi_deg(:))));
            else
                dphi_rad = 2 * pi;
            end
            sin_theta = sin(theta_rad);
            p_total   = ManualWeightsTuner.sph_power(af_copol, sin_theta, dtheta_rad, dphi_rad) + ...
                        ManualWeightsTuner.sph_power(af_cross,  sin_theta, dtheta_rad, dphi_rad);

            switch obj.active_polarization
                case ManualWeightsTuner.POLARIZATION_CROSS
                    power_linear_grid    = abs(af_cross) .^ 2;
                    metrics_array_factor = af_cross;
                case ManualWeightsTuner.POLARIZATION_TOTAL
                    power_linear_grid    = abs(af_copol) .^ 2 + abs(af_cross) .^ 2;
                    metrics_array_factor = sqrt(power_linear_grid);
                otherwise
                    power_linear_grid    = abs(af_copol) .^ 2;
                    metrics_array_factor = af_copol;
            end

            obj.dbi_grid = compute_directivity_dbi_grid( ...
                metrics_array_factor, obj.theta_deg, obj.phi_deg, p_total);

            if strcmp(obj.display_mode, ManualWeightsTuner.DISPLAY_ABSOLUTE)
                grid_display = obj.dbi_grid;
                clim_min = obj.dbi_min_field.Value;
                clim_max = obj.dbi_max_field.Value;
                if clim_max <= clim_min
                    clim_min = ManualWeightsTuner.DEFAULT_DBI_MIN;
                    clim_max = ManualWeightsTuner.DEFAULT_DBI_MAX;
                end
                obj.cbar.Ticks = linspace(clim_min, clim_max, 5);
            else
                power_db_grid = 10 * log10(max(power_linear_grid, ManualWeightsTuner.LOG10_EPSILON));
                grid_display  = max(power_db_grid - max(power_db_grid(:)), ...
                    -ManualWeightsTuner.DYNAMIC_RANGE_DB);
                clim_min = -ManualWeightsTuner.DYNAMIC_RANGE_DB;
                clim_max =  0.0;
            end

            obj.update_heatmap(grid_display, clim_min, clim_max);

            directives = obj.get_directives_from_data();
            obj.update_directive_overlays(directives);

            metrics = evaluate_metrics( ...
                obj.element_patterns_copol, obj.theta_deg, obj.phi_deg, ...
                obj.weights_complex, directives, [], metrics_array_factor, p_total);
            obj.update_metrics_labels(metrics);
        end

        function update_heatmap(obj, grid, clim_min, clim_max)
            obj.last_display_grid  = grid;
            obj.heatmap_surf.CData = grid;
            clim(obj.ax, [clim_min, clim_max]);
            drawnow limitrate;
        end

        function update_directive_overlays(obj, directives)
        % Remove old overlay artists; redraw contours + centre markers.
            ax_children = obj.ax.Children;
            to_delete   = ax_children(ax_children ~= obj.heatmap_surf & ...
                                      ax_children ~= obj.hover_text_h);
            delete(to_delete);

            if isempty(directives), return; end

            phys_masks = build_directive_physical_masks( ...
                obj.theta_deg, obj.phi_deg, directives);

            hold(obj.ax, 'on');
            for k = 1:numel(directives)
                d = directives{k};
                if strcmp(d.type, 'peak'), color = [0 0.55 0]; else, color = [0.85 0 0]; end

                mask = phys_masks{k};
                if any(mask(:))
                    contour(obj.ax, obj.phi_deg, obj.theta_deg, double(mask), ...
                        [0.5 0.5], 'Color', color, 'LineWidth', 2);
                end

                phi_d = 0.0;
                if isfield(d, 'phi'), phi_d = d.phi; end
                s = scatter(obj.ax, phi_d, d.theta, 120, color, 'Marker', '+');
                s.LineWidth = 2;
            end
            hold(obj.ax, 'off');
        end

        function update_metrics_labels(obj, metrics)
        % Refresh fixed metric labels and inline directive result labels.
        %
        % Uses char() for Greek letters because uilabel renders plain text
        % (not TeX/LaTeX), unlike uiaxes labels.
            th = char(952);   % θ
            ph = char(966);   % φ
            deg_s = char(176); % °

            obj.label_cost.Text = sprintf('%.4e', metrics.total_cost);
            obj.label_peak.Text = sprintf('%+.2f dBi', metrics.global_peak_dbi);
            obj.label_peak_angle.Text = sprintf('%s=%.1f%s  %s=%.1f%s', ...
                th, metrics.global_peak_theta_deg, deg_s, ...
                ph, metrics.global_peak_phi_deg,   deg_s);
            obj.label_hpbw.Text = sprintf('%s:%.1f%s  %s:%.1f%s', ...
                th, metrics.hpbw_theta_deg, deg_s, ...
                ph, metrics.hpbw_phi_deg,   deg_s);

            dm = metrics.directive_metrics;
            for k = 1:numel(obj.directive_widget_rows)
                lbl = obj.directive_widget_rows{k}.metric_lbl;
                if k > numel(dm)
                    lbl.Text      = char(8212);
                    lbl.FontColor = [0.4 0.4 0.4];
                    continue;
                end
                d_met = dm(k);
                if strcmp(d_met.type, 'peak')
                    lbl.Text      = sprintf('%+.2f dBi', d_met.gain_dbi);
                    lbl.FontColor = [0 0.55 0];
                else
                    lbl.Text      = sprintf('%+.2f dBi (%+.1f dBr)', ...
                        d_met.gain_dbi, d_met.null_depth_db);
                    lbl.FontColor = [0.85 0 0];
                end
            end
        end

    end % recompute methods

    % ══════════════════════════════════════════════════════════════
    % PRIVATE: HELPERS
    % ══════════════════════════════════════════════════════════════
    methods (Access = private)

        function directives = get_directives_from_data(obj)
            directives = {};
            for k = 1:numel(obj.directive_data)
                d = obj.directive_data{k};
                if ~ismember(d.type, {'peak', 'null'}), continue; end
                ds.type        = d.type;
                ds.theta       = d.theta;
                ds.phi         = d.phi;
                ds.theta_width = d.theta_width;
                ds.phi_width   = d.phi_width;
                ds.weight      = d.weight;
                directives{end + 1} = ds; %#ok<AGROW>
            end
        end

        function sync_entries_to_weights(obj)
        % Sync all editfields and sliders to match weights_complex.
            obj.syncing_weight_display = true;
            for n = 1:obj.n_elements
                amplitude_linear = abs(obj.weights_complex(n));
                phase_deg        = angle(obj.weights_complex(n)) * 180 / pi;
                obj.amp_fields{n}.Value    = amplitude_linear;
                obj.phase_fields{n}.Value  = phase_deg;
                obj.amp_sliders{n}.Value   = min(amplitude_linear, ManualWeightsTuner.AMP_SLIDER_MAX);
                obj.phase_sliders{n}.Value = max(ManualWeightsTuner.PHASE_SLIDER_MIN_DEG, ...
                    min(ManualWeightsTuner.PHASE_SLIDER_MAX_DEG, phase_deg));
            end
            obj.syncing_weight_display = false;
        end

        function weights = parse_weights_csv(obj, filepath)
        % Parse a weights CSV and return a complex weight vector (N_elements x 1).
            T = readtable(filepath, 'Delimiter', ',');
            if ~ismember('amplitude', T.Properties.VariableNames)
                error('ManualWeightsTuner:BadCSV', 'CSV missing ''amplitude'' column.');
            end
            if ~ismember('phase_deg', T.Properties.VariableNames)
                error('ManualWeightsTuner:BadCSV', 'CSV missing ''phase_deg'' column.');
            end
            if height(T) ~= obj.n_elements
                error('ManualWeightsTuner:CountMismatch', ...
                    ['CSV has %d row(s) but %d element(s) loaded. ' ...
                     'Ensure the CSV matches element_patterns_dir in config.'], ...
                    height(T), obj.n_elements);
            end
            weights = T.amplitude(:) .* exp(1i * deg2rad(T.phase_deg(:)));
        end

        function set_status(obj, message)
            obj.status_label.Text = sprintf('  %s', message);
        end

        function update_status_count(obj)
            obj.set_status(sprintf('%d elements loaded  %s  %d directives active', ...
                obj.n_elements, char(183), numel(obj.directive_data)));
        end

        function validate_config(obj)
            for key = {'element_patterns_dir', 'directives'}
                if ~isfield(obj.config, key{1})
                    error('ManualWeightsTuner:MissingKey', ...
                        'Required config key ''%s'' missing from config.yaml.', key{1});
                end
            end
        end

    end % helper methods

    % ══════════════════════════════════════════════════════════════
    % PRIVATE STATIC HELPERS
    % ══════════════════════════════════════════════════════════════
    methods (Access = private, Static)

        function p = sph_power(af, sin_theta, dtheta, dphi)
        % Spherical integral of |AF|^2 sin(θ) for total radiated power.
            p = sum(abs(af) .^ 2 .* sin_theta, 'all') * dtheta * dphi;
        end

        function b = dir_basename(p)
        % Return the last path component of a directory path string.
            p_clean = strtrim(p);
            while ~isempty(p_clean) && (p_clean(end) == '/' || p_clean(end) == '\')
                p_clean = p_clean(1:end - 1);
            end
            parts = strsplit(p_clean, {filesep, '/', '\'});
            b = parts{end};
            if isempty(b) && numel(parts) > 1, b = parts{end - 1}; end
        end

    end % static methods

end % classdef


% ── File-level helper: thin separator panel in a toolbar grid ─────
function mk_sep(parent_grid, col_index)
    p = uipanel(parent_grid, 'BorderType', 'none', 'BackgroundColor', [0.6 0.6 0.6]);
    p.Layout.Column = col_index;
end
