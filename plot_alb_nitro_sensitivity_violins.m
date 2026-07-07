function summary = plot_alb_nitro_sensitivity_violins(cases)
%PLOT_ALB_NITRO_SENSITIVITY_VIOLINS
% Violin plots for Alb + nitroimidazole structural sensitivity chains.
%
% Required input:
%   cases(k).chain_file
%   cases(k).regimen
%   cases(k).label
%
% Layout:
%   top-left     nitroimidazole Emax
%   top-right    nitroimidazole concentration scale
%   bottom-left  albendazole Emax
%   bottom-right albendazole concentration scale
%
% Each violin shows posterior log-multiplier samples after burn-in.

%% ========================================================================
%  USER TOGGLES
%  ========================================================================

cfg = struct();

cfg.burn_frac = 0.5;
cfg.font_size = 14;
cfg.figure_position = [100 100 1400 850];

cfg.violin_width = 0.34;
cfg.face_alpha = 0.28;
cfg.edge_alpha = 0.65;
cfg.line_width = 1.4;
cfg.median_line_width = 2.4;
cfg.interval_line_width = 1.4;

cfg.show_95_interval = true;
cfg.show_50_interval = true;
cfg.show_median = true;
cfg.show_baseline = true;

cfg.y_lim = [];     % e.g. [-2.5 2.5], or [] for automatic

cfg.show_titles = true;
cfg.show_legend = true;
cfg.legend_location = 'eastoutside';

% Colours: one colour per regimen, easy to change here
C = lines(max(7, numel(cases)));

for k = 1:numel(cases)
    if ~isfield(cases(k), 'color') || isempty(cases(k).color)
        cases(k).color = C(k,:);
    end
end

%% ========================================================================
%  Parameters to extract
%  ========================================================================

param_info = struct([]);

param_info(1).name = 'log_nitro_Emax';
param_info(1).title = 'Nitroimidazole E_{max}';
param_info(1).row = 1;
param_info(1).col = 1;

param_info(2).name = 'log_nitro_MIC';
param_info(2).title = 'Nitroimidazole concentration scale';
param_info(2).row = 1;
param_info(2).col = 2;

param_info(3).name = 'log_ABZ_Emax';
param_info(3).title = 'Albendazole E_{max}';
param_info(3).row = 2;
param_info(3).col = 1;

param_info(4).name = 'log_ABZ_MIC';
param_info(4).title = 'Albendazole concentration scale';
param_info(4).row = 2;
param_info(4).col = 2;

n_cases = numel(cases);
n_params = numel(param_info);

%% ========================================================================
%  Load samples
%  ========================================================================

samples = cell(n_cases, n_params);
summary = table();

for k = 1:n_cases

    chain = local_load_chain(cases(k).chain_file);

    Theta = chain.theta;
    burn_idx = floor(cfg.burn_frac * size(Theta,1)) + 1;
    Theta = Theta(burn_idx:end,:);

    names = chain.param_names;
    if isstring(names)
        names = cellstr(names);
    end

    for p = 1:n_params
        idx = find(strcmp(names, param_info(p).name), 1);

        if isempty(idx)
            error('Parameter %s not found in %s.', ...
                param_info(p).name, cases(k).chain_file);
        end

        s = Theta(:,idx);
        samples{k,p} = s;

        row = table();
        row.regimen = string(cases(k).regimen);
        row.label = string(cases(k).label);
        row.parameter = string(param_info(p).name);
        row.n_samples = numel(s);
        row.mean_log = mean(s);
        row.median_log = median(s);
        row.q025_log = prctile(s, 2.5);
        row.q25_log = prctile(s, 25);
        row.q75_log = prctile(s, 75);
        row.q975_log = prctile(s, 97.5);
        row.median_multiplier = exp(median(s));

        summary = [summary; row]; %#ok<AGROW>
    end
end

disp(summary)

%% ========================================================================
%  Plot
%  ========================================================================

fig = figure('Color','w','Position',cfg.figure_position);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax_all = gobjects(n_params,1);

for p = 1:n_params

    ax = nexttile;
    ax_all(p) = ax;
    hold(ax, 'on')

    for k = 1:n_cases

        s = samples{k,p};

        local_violin(ax, ...
            k, ...
            s, ...
            cases(k).color, ...
            cfg.violin_width, ...
            cfg.face_alpha, ...
            cfg.edge_alpha, ...
            cfg.line_width);

        q025 = prctile(s, 2.5);
        q25  = prctile(s, 25);
        q50  = prctile(s, 50);
        q75  = prctile(s, 75);
        q975 = prctile(s, 97.5);

        if cfg.show_95_interval
            plot(ax, [k k], [q025 q975], '-', ...
                'Color', cases(k).color, ...
                'LineWidth', cfg.interval_line_width, ...
                'HandleVisibility','off');
        end

        if cfg.show_50_interval
            plot(ax, [k k], [q25 q75], '-', ...
                'Color', cases(k).color, ...
                'LineWidth', cfg.interval_line_width + 1.2, ...
                'HandleVisibility','off');
        end

        if cfg.show_median
            plot(ax, [k-0.17 k+0.17], [q50 q50], '-', ...
                'Color', cases(k).color, ...
                'LineWidth', cfg.median_line_width, ...
                'HandleVisibility','off');
        end
    end

    if cfg.show_baseline
        yline(ax, 0, ':', 'Baseline', ...
            'Color',[0.25 0.25 0.25], ...
            'LineWidth',1.1, ...
            'HandleVisibility','off');
    end

    xlim(ax, [0.4 n_cases+0.6])

    set(ax, ...
        'XTick', 1:n_cases, ...
        'XTickLabel', {cases.label}, ...
        'FontSize', cfg.font_size);

    xtickangle(ax, 25)

    ylabel(ax, 'log-multiplier')
    grid(ax, 'on')
    box(ax, 'on')

    if cfg.show_titles
        title(ax, param_info(p).title, ...
            'FontSize', cfg.font_size + 1, ...
            'FontWeight','bold', ...
            'Interpreter','tex');
    end

    if ~isempty(cfg.y_lim)
        ylim(ax, cfg.y_lim)
    end
end

% Link y-limits by row, not globally
linkaxes(ax_all([1 2]), 'y')
linkaxes(ax_all([3 4]), 'y')

% Legend: dummy handles, one per regimen
if cfg.show_legend
    ax_leg = ax_all(2);
    h = gobjects(n_cases,1);

    for k = 1:n_cases
        h(k) = plot(ax_leg, NaN, NaN, '-', ...
            'Color', cases(k).color, ...
            'LineWidth', 3, ...
            'DisplayName', cases(k).label);
    end

    lgd = legend(ax_leg, h, 'Location', cfg.legend_location);
    lgd.FontSize = cfg.font_size;
end

% No global title by default; captions/slides can handle it
% sgtitle(tl, 'Albendazole + nitroimidazole structural sensitivity', ...
%     'FontSize', cfg.font_size + 2, 'FontWeight','bold');

end

%% ========================================================================
function chain = local_load_chain(chain_file)

S = load(chain_file);

if isfield(S, 'chain')
    chain = S.chain;
    return
end

fn = fieldnames(S);
for i = 1:numel(fn)
    candidate = S.(fn{i});
    if isstruct(candidate) && isfield(candidate, 'theta') && isfield(candidate, 'param_names')
        chain = candidate;
        return
    end
end

error('No chain struct found in %s.', chain_file);

end

%% ========================================================================
function local_violin(ax, x0, s, color, max_width, face_alpha, edge_alpha, line_width)

s = s(:);
s = s(isfinite(s));

if numel(s) < 5
    return
end

% KDE support
lo = prctile(s, 0.5);
hi = prctile(s, 99.5);
pad = 0.15 * max(1e-6, hi - lo);

ygrid = linspace(lo - pad, hi + pad, 250);

try
    f = ksdensity(s, ygrid);
catch
    % Fallback if ksdensity has issues
    [f, ygrid] = ksdensity(s);
end

if max(f) > 0
    f = f ./ max(f) * max_width;
end

x_left  = x0 - f;
x_right = x0 + f;

x_patch = [x_left fliplr(x_right)];
y_patch = [ygrid fliplr(ygrid)];

p = patch(ax, x_patch, y_patch, color, ...
    'FaceAlpha', face_alpha, ...
    'EdgeColor', color, ...
    'LineWidth', line_width, ...
    'HandleVisibility','off');

% Edge alpha not directly supported in older MATLAB patch edge;
% keep regular edge colour for compatibility.
if isprop(p, 'EdgeAlpha')
    p.EdgeAlpha = edge_alpha;
end

end

%{
C = lines(7);

cases_sens = struct([]);

% Canonical/reference first
cases_sens(1).chain_file = 'mcmc_outputs/pooled_MTZ14D_ABZ7D_4chains_currentSimulator.mat';
cases_sens(1).regimen = 'MTZ14D_ABZ7D';
cases_sens(1).scenario = 'refractory';
cases_sens(1).label = 'MTZ14D+ABZ7D';
cases_sens(1).color = C(1,:);

% Sensitivity one-chain runs
cases_sens(2).chain_file = 'mcmc_sensitivity_outputs/pooled_alb_nitro_mtz5_abz5_onechain_mcmc142_vp1001_currentSimulator.mat';
cases_sens(2).regimen = 'MTZ5D_ABZ5D';
cases_sens(2).scenario = 'refractory';
cases_sens(2).label = 'MTZ5D+ABZ5D';
cases_sens(2).color = C(2,:);

cases_sens(3).chain_file = 'mcmc_sensitivity_outputs/pooled_alb_nitro_mtz7_abz7_onechain_mcmc143_vp1001_currentSimulator.mat';
cases_sens(3).regimen = 'MTZ7D_ABZ7D';
cases_sens(3).scenario = 'refractory';
cases_sens(3).label = 'MTZ7D+ABZ7D';
cases_sens(3).color = C(3,:);

cases_sens(4).chain_file = 'mcmc_sensitivity_outputs/pooled_alb_nitro_tdz2g_abz5_onechain_mcmc144_vp1001_currentSimulator.mat';
cases_sens(4).regimen = 'TDZ2G_ABZ5D';
cases_sens(4).scenario = 'refractory';
cases_sens(4).label = 'TDZ2G+ABZ5D';
cases_sens(4).color = C(4,:);

cases_sens(5).chain_file = 'mcmc_sensitivity_outputs/pooled_alb_nitro_tdz2g_abz7_onechain_mcmc145_vp1001_currentSimulator.mat';
cases_sens(5).regimen = 'TDZ2G_ABZ7D';
cases_sens(5).scenario = 'refractory';
cases_sens(5).label = 'TDZ2G+ABZ7D';
cases_sens(5).color = C(5,:);

summary_violin = plot_alb_nitro_sensitivity_violins(cases_sens);
%}