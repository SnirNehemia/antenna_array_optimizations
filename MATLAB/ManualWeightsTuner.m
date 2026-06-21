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
        DEFAULT_DBI_MIN      = -30.0
        DEFAULT_DBI_MAX      =  10.0
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
        element_pattern_stacks
        polarization_options
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

        % Optimization tab widgets
        opt_max_iter_field
        opt_cost_tol_field
        opt_grad_tol_field
        opt_n_restarts_field
        opt_amp_min_field
        opt_amp_max_field
        opt_amp_unbounded_check
        opt_phase_only_check
        opt_amplitude_only_check
        opt_use_uniform_check
        opt_use_single_check
        opt_run_btn
        opt_status_label
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

            component_names = sort(fieldnames(patterns(1).components));
            obj.element_pattern_stacks = struct();
            for k = 1:numel(component_names)
                obj.element_pattern_stacks.(component_names{k}) = stack_component(patterns, component_names{k});
            end
            obj.polarization_options = [component_names(:); {ManualWeightsTuner.POLARIZATION_TOTAL}];

            obj.weights_complex        = ones(obj.n_elements, 1, 'like', 1 + 0i);
            obj.active_polarization    = ManualWeightsTuner.POLARIZATION_TOTAL;
            obj.display_mode           = ManualWeightsTuner.DISPLAY_ABSOLUTE;
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

            % Calls on_display_mode_change (not recompute_and_redraw directly) so the
            % colorbar label/ticks are correctly initialized for the default
            % absolute-display mode before the first render.
            obj.on_display_mode_change();
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

            % Top-level 2-row grid: [toolbar (fixed); content (fills)].
            top_grid = uigridlayout(obj.fig, [2, 1]);
            top_grid.RowHeight   = {96, '1x'};
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

            % Right column: 2-row grid directly in content_grid.
            right_grid = uigridlayout(content_grid, [2, 1]);
            right_grid.Layout.Row    = 1;
            right_grid.Layout.Column = 2;
            right_grid.RowHeight     = {ManualWeightsTuner.WEIGHTS_PANEL_HEIGHT + 44, '1x'};
            right_grid.ColumnWidth   = {'1x'};
            right_grid.Padding       = [0 0 0 0];
            right_grid.RowSpacing    = 4;

            obj.build_weights_panel(right_grid);

            % Tab group: "Directives" | "Optimization" (run optimizer) | "Metrics"
            tabgroup = uitabgroup(right_grid);
            tabgroup.Layout.Row    = 2;
            tabgroup.Layout.Column = 1;

            tab_directives = uitab(tabgroup, 'Title', 'Directives');
            obj.build_directives_panel(tab_directives);

            tab_optimization = uitab(tabgroup, 'Title', 'Optimization');
            obj.build_optimization_panel(tab_optimization);

            tab_metrics = uitab(tabgroup, 'Title', 'Metrics');
            obj.build_metrics_panel(tab_metrics);
        end

        % ── BUILD TOOLBAR ─────────────────────────────────────────

        function build_toolbar(obj, parent)
        % Build a 2-row toolbar:
        %   Row 1: info block | Load Config | Load Data Folder | [sep] |
        %          Display combo | dBi min/max | [sep] | Polarization combo
        %   Row 2: Uniform Weights | Load Weights CSV
        %
        % Labels intentionally have NO explicit BackgroundColor so MATLAB's
        % theme (including dark mode) renders them correctly.
            outer = uigridlayout(parent, [2, 1]);
            outer.RowHeight   = {'1x', '1x'};
            outer.ColumnWidth = {'1x'};
            outer.Padding     = [6 4 6 4];
            outer.RowSpacing  = 2;

            % ══ Row 1 ═══════════════════════════════════════════════
            tg = uigridlayout(outer, [1, 13]);
            tg.Layout.Row    = 1;
            tg.Layout.Column = 1;
            tg.ColumnWidth   = {'1x', 100, 110, 10, 60, 90, 58, 64, 58, 64, 10, 80, 110};
            tg.Padding       = [0 0 0 0];
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

            % ── Load Data Folder button ────────────────────────────
            btn_data = uibutton(tg, 'Text', 'Load Data Folder...', ...
                'FontSize', 11, ...
                'ButtonPushedFcn', @(~,~) obj.on_load_data_folder());
            btn_data.Layout.Column = 3;

            % separator (thin panel)
            mk_sep(tg, 4);

            % ── Display mode ──────────────────────────────────────
            lbl_disp = uilabel(tg, 'Text', 'Display:', 'FontSize', 11, ...
                'HorizontalAlignment', 'right');
            lbl_disp.Layout.Column = 5;

            obj.disp_dropdown = uidropdown(tg, ...
                'Items',           {ManualWeightsTuner.DISPLAY_RELATIVE, ManualWeightsTuner.DISPLAY_ABSOLUTE}, ...
                'Value',           ManualWeightsTuner.DISPLAY_ABSOLUTE, ...
                'FontSize',        11, ...
                'ValueChangedFcn', @(~,~) obj.on_display_mode_change());
            obj.disp_dropdown.Layout.Column = 6;

            lbl_min = uilabel(tg, 'Text', 'min (dBi):', 'FontSize', 11, ...
                'HorizontalAlignment', 'right');
            lbl_min.Layout.Column = 7;

            obj.dbi_min_field = uieditfield(tg, 'numeric', ...
                'Value',           ManualWeightsTuner.DEFAULT_DBI_MIN, ...
                'FontSize',        11, ...
                'ValueChangedFcn', @(~,~) obj.on_display_mode_change());
            obj.dbi_min_field.Layout.Column = 8;

            lbl_max = uilabel(tg, 'Text', 'max (dBi):', 'FontSize', 11, ...
                'HorizontalAlignment', 'right');
            lbl_max.Layout.Column = 9;

            obj.dbi_max_field = uieditfield(tg, 'numeric', ...
                'Value',           ManualWeightsTuner.DEFAULT_DBI_MAX, ...
                'FontSize',        11, ...
                'ValueChangedFcn', @(~,~) obj.on_display_mode_change());
            obj.dbi_max_field.Layout.Column = 10;

            mk_sep(tg, 11);

            % ── Polarization ──────────────────────────────────────
            lbl_pol = uilabel(tg, 'Text', 'Polarization:', 'FontSize', 11, ...
                'HorizontalAlignment', 'right');
            lbl_pol.Layout.Column = 12;

            obj.pol_dropdown = uidropdown(tg, ...
                'Items',           obj.polarization_options, ...
                'Value',           obj.active_polarization, ...
                'FontSize',        11, ...
                'ValueChangedFcn', @(~,~) obj.on_polarization_change());
            obj.pol_dropdown.Layout.Column = 13;

            % ══ Row 2 ═══════════════════════════════════════════════
            bg = uigridlayout(outer, [1, 3]);
            bg.Layout.Row    = 2;
            bg.Layout.Column = 1;
            bg.ColumnWidth   = {'1x', 140, 140};
            bg.Padding       = [0 0 0 0];
            bg.ColumnSpacing = 6;

            btn_uni = uibutton(bg, 'Text', 'Uniform Weights', ...
                'FontSize', 11, ...
                'ButtonPushedFcn', @(~,~) obj.on_uniform());
            btn_uni.Layout.Column = 2;

            btn_csv = uibutton(bg, 'Text', 'Load Weights CSV', ...
                'FontSize', 11, ...
                'ButtonPushedFcn', @(~,~) obj.on_load_csv());
            btn_csv.Layout.Column = 3;
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

            % uigridlayout children are always resized to fill their parent
            % (no natural/overflow size), so a uigridlayout cannot trigger
            % scrolling inside a Scrollable uipanel. Instead, lay out the
            % rows with absolute pixel positions on a fixed-height content
            % panel that is taller than the visible scroll area.
            n_rows = obj.n_elements + 1;   % +1 for header
            ROW_H   = 26;
            ROW_GAP = 2;
            PAD     = 4;
            W = ManualWeightsTuner.RIGHT_PANEL_WIDTH - 24;
            H = PAD * 2 + n_rows * ROW_H + (n_rows - 1) * ROW_GAP;

            content = uipanel(outer, 'BorderType', 'none', ...
                'Position', [1 1 W H], 'AutoResizeChildren', 'off');

            % Column layout (fixed pixel widths, 3px spacing).
            w_idx = 40; w_amp_ef = 62; w_phase_ef = 66; w_solo = 46;
            w_slider = (W - 2 * PAD - w_idx - w_amp_ef - w_phase_ef - w_solo - 5 * 3) / 2;
            x_idx      = PAD;
            x_amp_ef   = x_idx + w_idx + 3;
            x_amp_sl   = x_amp_ef + w_amp_ef + 3;
            x_phase_ef = x_amp_sl + w_slider + 3;
            x_phase_sl = x_phase_ef + w_phase_ef + 3;
            x_solo     = x_phase_sl + w_slider + 3;

            row_y = @(r) H - PAD - r * ROW_H - (r - 1) * ROW_GAP;

            % Header row.
            headers = {'Elem', 'Amplitude', '', 'Phase (°)', '', ''};
            xs = [x_idx, x_amp_ef, x_amp_sl, x_phase_ef, x_phase_sl, x_solo];
            ws = [w_idx, w_amp_ef, w_slider, w_phase_ef, w_slider, w_solo];
            y = row_y(1);
            for c = 1:6
                uilabel(content, 'Text', headers{c}, ...
                    'FontWeight', 'bold', 'FontSize', 10, ...
                    'HorizontalAlignment', 'center', ...
                    'Position', [xs(c) y ws(c) ROW_H]);
            end

            for n = 1:obj.n_elements
                y = row_y(n + 1);

                uilabel(content, 'Text', sprintf('%3d', n - 1), ...
                    'FontSize', 10, 'HorizontalAlignment', 'right', ...
                    'Position', [x_idx y w_idx ROW_H]);

                amp_ef = uieditfield(content, 'numeric', ...
                    'Value',           ManualWeightsTuner.INITIAL_AMPLITUDE, ...
                    'Limits',          [0 Inf], ...
                    'FontSize',        10, ...
                    'Position',        [x_amp_ef y w_amp_ef ROW_H], ...
                    'ValueChangedFcn', @(src,~) obj.on_amp_entry_changed(n, src.Value));

                amp_sl = uislider(content, ...
                    'Limits',           [ManualWeightsTuner.AMP_SLIDER_MIN, ManualWeightsTuner.AMP_SLIDER_MAX], ...
                    'Value',            ManualWeightsTuner.INITIAL_AMPLITUDE, ...
                    'MajorTicks',       [], 'MinorTicks', [], ...
                    'Position',         [x_amp_sl y+ROW_H/2 w_slider 3], ...
                    'ValueChangedFcn',  @(src,~) obj.on_amp_slider_changed(n, src.Value), ...
                    'ValueChangingFcn', @(src,evt) obj.on_amp_slider_changing(n, evt.Value));

                phase_ef = uieditfield(content, 'numeric', ...
                    'Value',           ManualWeightsTuner.INITIAL_PHASE_DEG, ...
                    'FontSize',        10, ...
                    'Position',        [x_phase_ef y w_phase_ef ROW_H], ...
                    'ValueChangedFcn', @(src,~) obj.on_phase_entry_changed(n, src.Value));

                phase_sl = uislider(content, ...
                    'Limits',           [ManualWeightsTuner.PHASE_SLIDER_MIN_DEG, ManualWeightsTuner.PHASE_SLIDER_MAX_DEG], ...
                    'Value',            ManualWeightsTuner.INITIAL_PHASE_DEG, ...
                    'MajorTicks',       [], 'MinorTicks', [], ...
                    'Position',         [x_phase_sl y+ROW_H/2 w_slider 3], ...
                    'ValueChangedFcn',  @(src,~) obj.on_phase_slider_changed(n, src.Value), ...
                    'ValueChangingFcn', @(src,evt) obj.on_phase_slider_changing(n, evt.Value));

                uibutton(content, 'Text', 'Solo', ...
                    'FontSize',        10, ...
                    'Position',        [x_solo y w_solo ROW_H], ...
                    'ButtonPushedFcn', @(~,~) obj.on_solo_element(n));

                obj.amp_fields{n}    = amp_ef;
                obj.amp_sliders{n}   = amp_sl;
                obj.phase_fields{n}  = phase_ef;
                obj.phase_sliders{n} = phase_sl;
            end
        end

        % ── BUILD DIRECTIVES PANEL ────────────────────────────────

        function build_directives_panel(obj, parent)
        % Tab content with dynamic add/remove directive rows.
            outer_grid = uigridlayout(parent, [3, 1]);
            outer_grid.RowHeight   = {24, '1x', 32};
            outer_grid.ColumnWidth = {'1x'};
            outer_grid.Padding     = [2 2 2 2];
            outer_grid.RowSpacing  = 2;

            % Column header row.
            th = char(952); ph = char(966);
            hdr_g = uigridlayout(outer_grid, [1, 9]);
            hdr_g.Layout.Row    = 1;
            hdr_g.Layout.Column = 1;
            hdr_g.ColumnWidth   = {52, 42, 42, 42, 42, 34, 64, 24, '1x'};
            hdr_g.Padding       = [0 0 0 0];
            hdr_g.ColumnSpacing = 2;
            hdr_names = {'Type', [th,'(°)'], [ph,'(°)'], [th,'W'], [ph,'W'], 'Wt', 'Agg', '', 'Result'};
            for c = 1:9
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

        % ── BUILD OPTIMIZATION PANEL ──────────────────────────────

        function build_optimization_panel(obj, parent)
        % Tab content: scrollable, editable optimizer settings (from
        % config.optimizer), each with a short description, plus a
        % "Run Optimization" button pinned at the bottom that runs
        % run_optimization() with the current directives/settings and loads
        % the resulting weights.
            opt_cfg = struct();
            if isfield(obj.config, 'optimizer')
                opt_cfg = obj.config.optimizer;
            end

            % Outer: scrollable settings area (fills) + pinned run bar (fixed).
            outer = uigridlayout(parent, [2, 1]);
            outer.RowHeight   = {'1x', 44};
            outer.ColumnWidth = {'1x'};
            outer.Padding     = [0 0 0 0];
            outer.RowSpacing  = 4;

            scroll_panel = uipanel(outer, 'BorderType', 'none', 'Scrollable', 'on');
            scroll_panel.Layout.Row    = 1;
            scroll_panel.Layout.Column = 1;

            % uigridlayout children are always resized to fill their parent
            % (no natural/overflow size), so a uigridlayout cannot trigger
            % scrolling inside a Scrollable uipanel. Instead, lay out the
            % settings with absolute pixel positions on a fixed-height
            % content panel that is taller than the visible scroll area.
            W        = ManualWeightsTuner.RIGHT_PANEL_WIDTH - 24;
            ROW_H    = 22;
            DESC_H   = 34;
            ROW_GAP  = 4;
            SEC_GAP  = 8;
            TOP_PAD  = 8;
            n_settings = 10;
            H = TOP_PAD + n_settings * (ROW_H + ROW_GAP + DESC_H + SEC_GAP);

            content = uipanel(scroll_panel, 'BorderType', 'none', ...
                'Position', [1 1 W H], 'AutoResizeChildren', 'off');

            x_label  = 8;   w_label = 130;
            x_field1 = x_label + w_label + 8;  w_field = 80;
            x_field2 = x_field1 + w_field + 8;
            x_full   = 8;   w_full  = W - 16;

            y = H - TOP_PAD;

            % Max iterations
            y = y - ROW_H;
            mk_label_pos(content, x_label, y, w_label, ROW_H, 'Max iterations:');
            obj.opt_max_iter_field = uieditfield(content, 'numeric', ...
                'Position', [x_field1 y w_field ROW_H], ...
                'Value', cfg_get(opt_cfg, 'max_iterations', 250), ...
                'Limits', [1 Inf], 'RoundFractionalValues', 'on', 'FontSize', 10);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Maximum L-BFGS-B iterations per restart. Stops earlier if cost_tolerance is met.');
            y = y - SEC_GAP;

            % Cost tolerance
            y = y - ROW_H;
            mk_label_pos(content, x_label, y, w_label, ROW_H, 'Cost tolerance:');
            obj.opt_cost_tol_field = uieditfield(content, 'numeric', ...
                'Position', [x_field1 y w_field ROW_H], ...
                'Value', cfg_get(opt_cfg, 'cost_tolerance', 1e-6), ...
                'Limits', [0 Inf], 'FontSize', 10);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Relative cost-improvement stopping threshold. Smaller = finer convergence, slower.');
            y = y - SEC_GAP;

            % Gradient tolerance
            y = y - ROW_H;
            mk_label_pos(content, x_label, y, w_label, ROW_H, 'Gradient tolerance:');
            obj.opt_grad_tol_field = uieditfield(content, 'numeric', ...
                'Position', [x_field1 y w_field ROW_H], ...
                'Value', cfg_get(opt_cfg, 'gradient_tolerance', 1e-5), ...
                'Limits', [0 Inf], 'FontSize', 10);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Gradient-norm stopping threshold. Whichever of this or cost_tolerance fires first wins.');
            y = y - SEC_GAP;

            % N restarts
            y = y - ROW_H;
            mk_label_pos(content, x_label, y, w_label, ROW_H, 'N restarts:');
            obj.opt_n_restarts_field = uieditfield(content, 'numeric', ...
                'Position', [x_field1 y w_field ROW_H], ...
                'Value', cfg_get(opt_cfg, 'n_restarts', 2), ...
                'Limits', [0 Inf], 'RoundFractionalValues', 'on', 'FontSize', 10);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Number of independent multi-start runs (plus one per element if single-element init is enabled).');
            y = y - SEC_GAP;

            % Amplitude bounds [min, max]
            y = y - ROW_H;
            mk_label_pos(content, x_label, y, w_label, ROW_H, 'Amplitude bounds:');
            amp_bounds = cfg_get(opt_cfg, 'amplitude_bounds', [0.0 1.0]);
            is_unbounded = isempty(amp_bounds);
            if is_unbounded, amp_bounds = [0.0 1.0]; end
            obj.opt_amp_min_field = uieditfield(content, 'numeric', ...
                'Position', [x_field1 y w_field ROW_H], ...
                'Value', amp_bounds(1), 'Limits', [0 Inf], 'FontSize', 10, ...
                'Enable', ~is_unbounded);
            obj.opt_amp_max_field = uieditfield(content, 'numeric', ...
                'Position', [x_field2 y w_field ROW_H], ...
                'Value', amp_bounds(2), 'Limits', [0 Inf], 'FontSize', 10, ...
                'Enable', ~is_unbounded);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Per-element amplitude [min, max] range. Ignored when phase_only is enabled.');
            y = y - SEC_GAP;

            % Unbounded checkbox
            y = y - ROW_H;
            obj.opt_amp_unbounded_check = uicheckbox(content, 'Text', 'Unbounded amplitude', ...
                'Position', [x_label y w_full ROW_H], ...
                'Value', is_unbounded, 'FontSize', 10, ...
                'ValueChangedFcn', @(src,~) obj.on_amp_unbounded_changed(src.Value));
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'When checked, amplitudes are left unbounded (amplitude_bounds: null).');
            y = y - SEC_GAP;

            % Phase only
            y = y - ROW_H;
            obj.opt_phase_only_check = uicheckbox(content, 'Text', 'Phase only', ...
                'Position', [x_label y w_full ROW_H], ...
                'Value', cfg_get(opt_cfg, 'phase_only', false), 'FontSize', 10);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Fix all amplitudes to 1.0 and optimize phase only. Overrides amplitude_bounds.');
            y = y - SEC_GAP;

            % Amplitude only
            y = y - ROW_H;
            obj.opt_amplitude_only_check = uicheckbox(content, 'Text', 'Amplitude only', ...
                'Position', [x_label y w_full ROW_H], ...
                'Value', cfg_get(opt_cfg, 'amplitude_only', false), 'FontSize', 10);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Fix all phases to 0 and optimize amplitudes only. Mutually exclusive with phase_only.');
            y = y - SEC_GAP;

            % Use uniform init
            y = y - ROW_H;
            obj.opt_use_uniform_check = uicheckbox(content, 'Text', 'Use uniform init', ...
                'Position', [x_label y w_full ROW_H], ...
                'Value', cfg_get(opt_cfg, 'use_uniform_init', true), 'FontSize', 10);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Include a deterministic uniform-weights run (amplitude=1, phase=0) as the first restart.');
            y = y - SEC_GAP;

            % Use single-element init
            y = y - ROW_H;
            obj.opt_use_single_check = uicheckbox(content, 'Text', 'Use single-element init', ...
                'Position', [x_label y w_full ROW_H], ...
                'Value', cfg_get(opt_cfg, 'use_single_element_init', true), 'FontSize', 10);
            y = y - ROW_GAP - DESC_H;
            mk_desc_pos(content, x_full, y, w_full, DESC_H, ...
                'Add one extra restart per element, starting with only that element active (adds N_elements restarts).');

            % ── Pinned run bar ─────────────────────────────────────
            run_grid = uigridlayout(outer, [1, 2]);
            run_grid.Layout.Row    = 2;
            run_grid.Layout.Column = 1;
            run_grid.ColumnWidth   = {150, '1x'};
            run_grid.Padding       = [8 4 8 4];
            run_grid.ColumnSpacing = 8;

            obj.opt_run_btn = uibutton(run_grid, 'Text', 'Run Optimization', ...
                'FontSize', 11, 'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~,~) obj.on_run_optimization());
            obj.opt_run_btn.Layout.Row    = 1;
            obj.opt_run_btn.Layout.Column = 1;

            obj.opt_status_label = uilabel(run_grid, 'Text', '', ...
                'FontSize', 10, 'FontColor', [0.45 0.45 0.45]);
            obj.opt_status_label.Layout.Row    = 1;
            obj.opt_status_label.Layout.Column = 2;
        end

        % ── BUILD METRICS PANEL ───────────────────────────────────

        function build_metrics_panel(obj, parent)
        % Four fixed metric rows.
            met_grid = uigridlayout(parent, [4, 2]);
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
                       'theta_width', 5.0, 'phi_width', 5.0, 'weight', 1.0, ...
                       'aggregation', 'mean');
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
                if isfield(d, 'aggregation') && ~isempty(d.aggregation)
                    D.aggregation = d.aggregation;
                end
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
            row_grid = uigridlayout(obj.directives_inner_grid, [1, 9]);
            row_grid.Layout.Row    = row_k;
            row_grid.Layout.Column = 1;
            row_grid.ColumnWidth   = {52, 42, 42, 42, 42, 34, 64, 24, '1x'};
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

            agg_dd = uidropdown(row_grid, ...
                'Items',           {'mean', 'max', 'min'}, ...
                'Value',           D.aggregation, ...
                'FontSize',        10, ...
                'ValueChangedFcn', @(src,~) obj.on_directive_field_changed(row_k, 'aggregation', src.Value));
            agg_dd.Layout.Column = 7;

            rm_btn = uibutton(row_grid, 'Text', char(215), ...
                'FontSize',        12, ...
                'ButtonPushedFcn', @(~,~) obj.remove_directive_row(row_k));
            rm_btn.Layout.Column = 8;

            metric_lbl = uilabel(row_grid, 'Text', char(8212), ...
                'FontSize', 10, 'FontName', 'Courier New');
            metric_lbl.Layout.Column = 9;

            rw.row_grid   = row_grid;
            rw.type_dd    = type_dd;
            rw.theta_ef   = efs{1};
            rw.phi_ef     = efs{2};
            rw.tw_ef      = efs{3};
            rw.pw_ef      = efs{4};
            rw.weight_ef  = efs{5};
            rw.agg_dd     = agg_dd;
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

        function on_load_data_folder(obj)
        % Pick a new element_patterns_dir, write a temp config with that
        % folder substituted, and reload via the next_config_path mechanism.
            patterns_dir = obj.config.element_patterns_dir;
            if ~isfolder(patterns_dir)
                repo_root    = fileparts(fileparts(mfilename('fullpath')));
                patterns_dir = fullfile(repo_root, patterns_dir);
            end
            start_dir = fileparts(patterns_dir);
            folder = uigetdir(start_dir, 'Select Element Patterns Data Folder');
            if isequal(folder, 0), return; end

            raw = fileread(obj.config_path_str);
            new_line = sprintf('element_patterns_dir: "%s"', strrep(folder, '\', '/'));
            rawlines = regexp(raw, '\r\n|\r|\n', 'split');
            for i = 1:numel(rawlines)
                if ~isempty(regexp(rawlines{i}, '^element_patterns_dir:', 'once'))
                    rawlines{i} = new_line;
                    break;
                end
            end
            raw = strjoin(rawlines, newline);

            tmp_path = fullfile(tempdir, sprintf('manual_tuner_config_%s.yaml', ...
                char(datetime('now', 'Format', 'yyyyMMdd_HHmmssSSS'))));
            fid = fopen(tmp_path, 'w');
            fwrite(fid, raw);
            fclose(fid);

            obj.next_config_path = tmp_path;
            delete(obj.fig);
        end

        function on_amp_unbounded_changed(obj, is_unbounded)
            obj.opt_amp_min_field.Enable = ~is_unbounded;
            obj.opt_amp_max_field.Enable = ~is_unbounded;
        end

        function on_run_optimization(obj)
        % Build a temp config from the current directives + optimizer settings
        % and run the full run_optimization() pipeline, then load the
        % resulting weights into the tuner.
            directives = obj.get_directives_from_data();
            if isempty(directives)
                uialert(obj.fig, 'Add at least one directive before running the optimizer.', ...
                    'No Directives');
                return;
            end

            optimizer_cfg = struct();
            optimizer_cfg.max_iterations      = obj.opt_max_iter_field.Value;
            optimizer_cfg.cost_tolerance      = obj.opt_cost_tol_field.Value;
            optimizer_cfg.gradient_tolerance  = obj.opt_grad_tol_field.Value;
            optimizer_cfg.n_restarts          = obj.opt_n_restarts_field.Value;
            if obj.opt_amp_unbounded_check.Value
                optimizer_cfg.amplitude_bounds = [];
            else
                optimizer_cfg.amplitude_bounds = ...
                    [obj.opt_amp_min_field.Value, obj.opt_amp_max_field.Value];
            end
            optimizer_cfg.phase_only              = obj.opt_phase_only_check.Value;
            optimizer_cfg.amplitude_only          = obj.opt_amplitude_only_check.Value;
            optimizer_cfg.use_uniform_init        = obj.opt_use_uniform_check.Value;
            optimizer_cfg.use_single_element_init = obj.opt_use_single_check.Value;

            raw = fileread(obj.config_path_str);
            raw = ManualWeightsTuner.replace_yaml_section(raw, 'directives', ...
                ManualWeightsTuner.format_directives_yaml(directives));
            raw = ManualWeightsTuner.replace_yaml_section(raw, 'optimizer', ...
                ManualWeightsTuner.format_optimizer_yaml(optimizer_cfg));

            % Run with the polarization currently shown in the toolbar,
            % regardless of what config.yaml originally specified.
            pol_line = sprintf('polarization: "%s"', obj.active_polarization);
            rawlines = regexp(raw, '\r\n|\r|\n', 'split');
            found = false;
            for i = 1:numel(rawlines)
                if ~isempty(regexp(rawlines{i}, '^polarization:', 'once'))
                    rawlines{i} = pol_line;
                    found = true;
                    break;
                end
            end
            if ~found
                rawlines{end+1} = pol_line;
            end
            raw = strjoin(rawlines, newline);

            tmp_path = fullfile(tempdir, sprintf('manual_tuner_run_config_%s.yaml', ...
                char(datetime('now', 'Format', 'yyyyMMdd_HHmmssSSS'))));
            fid = fopen(tmp_path, 'w');
            fwrite(fid, raw);
            fclose(fid);

            obj.opt_run_btn.Enable = 'off';
            obj.opt_status_label.Text = 'Running optimization... (see MATLAB console for progress)';
            d = uiprogressdlg(obj.fig, 'Title', 'Optimizing', ...
                'Message', 'Running multi-start optimization... see MATLAB console for progress.', ...
                'Indeterminate', 'on');
            drawnow;

            try
                output_dir = run_optimization(tmp_path);
                weights_csv = fullfile(output_dir, 'weights.csv');
                obj.weights_complex = obj.parse_weights_csv(weights_csv);
                obj.sync_entries_to_weights();
                obj.recompute_and_redraw();
                obj.opt_status_label.Text = sprintf('Done. Results: %s', output_dir);
                obj.set_status(sprintf('Optimization complete %s results saved to: %s', char(183), output_dir));
            catch exc
                obj.opt_status_label.Text = 'Optimization failed.';
                close(d);
                uialert(obj.fig, exc.message, 'Optimization Error');
                obj.opt_run_btn.Enable = 'on';
                return;
            end

            close(d);
            obj.opt_run_btn.Enable = 'on';
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

            component_names = fieldnames(obj.element_pattern_stacks);
            af_components = struct();
            for k = 1:numel(component_names)
                af_components.(component_names{k}) = ...
                    compute_array_factor(w_norm, obj.element_pattern_stacks.(component_names{k}));
            end

            % Spherical total radiated power (CST partial-directivity normaliser).
            theta_rad = deg2rad(obj.theta_deg(:));
            if numel(theta_rad) > 1, dtheta_rad = mean(diff(theta_rad)); else, dtheta_rad = pi; end
            if numel(obj.phi_deg) > 1
                dphi_rad = deg2rad(mean(diff(obj.phi_deg(:))));
            else
                dphi_rad = 2 * pi;
            end
            sin_theta = sin(theta_rad);
            p_total = 0;
            for k = 1:numel(component_names)
                p_total = p_total + ManualWeightsTuner.sph_power( ...
                    af_components.(component_names{k}), sin_theta, dtheta_rad, dphi_rad);
            end

            if strcmp(obj.active_polarization, ManualWeightsTuner.POLARIZATION_TOTAL)
                power_linear_grid = zeros(size(obj.theta_deg, 1), numel(obj.phi_deg));
                for k = 1:numel(component_names)
                    power_linear_grid = power_linear_grid + abs(af_components.(component_names{k})) .^ 2;
                end
                metrics_array_factor = sqrt(power_linear_grid);
            else
                af = af_components.(obj.active_polarization);
                power_linear_grid    = abs(af) .^ 2;
                metrics_array_factor = af;
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

            ep_reference = obj.element_pattern_stacks.(component_names{1});
            metrics = evaluate_metrics( ...
                ep_reference, obj.theta_deg, obj.phi_deg, ...
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
                ds.aggregation = d.aggregation;
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

        function text = format_directives_yaml(directives)
        % Serialize a cell array of directive structs (as returned by
        % get_directives_from_data) into a YAML "directives:" block
        % compatible with read_config_yaml.
            lines = {'directives:'};
            for k = 1:numel(directives)
                d = directives{k};
                lines{end + 1} = sprintf('  - type: "%s"', d.type); %#ok<AGROW>
                lines{end + 1} = sprintf('    theta: %g', d.theta); %#ok<AGROW>
                lines{end + 1} = sprintf('    phi: %g', d.phi); %#ok<AGROW>
                lines{end + 1} = sprintf('    theta_width: %g', d.theta_width); %#ok<AGROW>
                lines{end + 1} = sprintf('    phi_width: %g', d.phi_width); %#ok<AGROW>
                lines{end + 1} = sprintf('    weight: %g', d.weight); %#ok<AGROW>
                lines{end + 1} = sprintf('    aggregation: "%s"', d.aggregation); %#ok<AGROW>
            end
            text = strjoin(lines, newline);
        end

        function text = format_optimizer_yaml(opt_cfg)
        % Serialize the optimizer settings struct into a YAML "optimizer:"
        % block compatible with read_config_yaml.
            lines = {'optimizer:'};
            lines{end + 1} = sprintf('  max_iterations: %d', round(opt_cfg.max_iterations));
            lines{end + 1} = sprintf('  cost_tolerance: %g', opt_cfg.cost_tolerance);
            lines{end + 1} = sprintf('  gradient_tolerance: %g', opt_cfg.gradient_tolerance);
            lines{end + 1} = sprintf('  n_restarts: %d', round(opt_cfg.n_restarts));
            if isempty(opt_cfg.amplitude_bounds)
                lines{end + 1} = '  amplitude_bounds: null';
            else
                lines{end + 1} = sprintf('  amplitude_bounds: [%g, %g]', ...
                    opt_cfg.amplitude_bounds(1), opt_cfg.amplitude_bounds(2));
            end
            lines{end + 1} = sprintf('  phase_only: %s', ManualWeightsTuner.bool_str(opt_cfg.phase_only));
            lines{end + 1} = sprintf('  amplitude_only: %s', ManualWeightsTuner.bool_str(opt_cfg.amplitude_only));
            lines{end + 1} = sprintf('  use_uniform_init: %s', ManualWeightsTuner.bool_str(opt_cfg.use_uniform_init));
            lines{end + 1} = sprintf('  use_single_element_init: %s', ManualWeightsTuner.bool_str(opt_cfg.use_single_element_init));
            text = strjoin(lines, newline);
        end

        function s = bool_str(tf)
            if tf, s = 'true'; else, s = 'false'; end
        end

        function raw = replace_yaml_section(raw, key, replacement_text)
        % Replace the top-level "key:" block (the key line through the line
        % before the next top-level "name:" key) in raw YAML text with
        % replacement_text. Used to splice updated directives/optimizer
        % blocks into a temp copy of config.yaml before running the optimizer.
            rawlines = regexp(raw, '\r\n|\r|\n', 'split');
            start_idx = find(~cellfun(@isempty, regexp(rawlines, ['^' key ':'], 'once')), 1);
            if isempty(start_idx)
                % Key not present: append the new block at the end.
                raw = [raw, newline, replacement_text, newline];
                return;
            end
            end_idx = numel(rawlines) + 1;
            for i = (start_idx + 1):numel(rawlines)
                if ~isempty(regexp(rawlines{i}, '^[A-Za-z_][A-Za-z0-9_]*:', 'once'))
                    end_idx = i;
                    break;
                end
            end
            new_lines = [rawlines(1:start_idx - 1), strsplit(replacement_text, newline), rawlines(end_idx:end)];
            raw = strjoin(new_lines, newline);
        end

    end % static methods

end % classdef


% ── File-level helper: thin separator panel in a toolbar grid ─────
function mk_sep(parent_grid, col_index)
    p = uipanel(parent_grid, 'BorderType', 'none', 'BackgroundColor', [0.6 0.6 0.6]);
    p.Layout.Column = col_index;
end

% ── File-level helper: right-aligned label at an absolute position ──
function lbl = mk_label_pos(parent, x, y, w, h, text)
    lbl = uilabel(parent, 'Text', text, 'FontSize', 10, ...
        'HorizontalAlignment', 'right', 'Position', [x y w h]);
end

% ── File-level helper: wrapped description label at an absolute position ──
function lbl = mk_desc_pos(parent, x, y, w, h, text)
    lbl = uilabel(parent, 'Text', text, 'FontSize', 9, ...
        'FontColor', [0.5 0.5 0.5], 'WordWrap', 'on', ...
        'VerticalAlignment', 'top', 'Position', [x y w h]);
end

% ── File-level helper: read optional config field with default ────
function val = cfg_get(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = default;
    end
end
