function summary = plot_qnc_prior_sensitivity_violins(cases)
%PLOT_QNC_PRIOR_SENSITIVITY_VIOLINS
% Prior sensitivity plot for QNC chains.
%
% INPUT:
%   cases(k).chain_file
%   cases(k).label
%
% Optional:
%   cases(k).color
%
% Assumes each chain contains:
%   log_QNC_Emax
%   log_QNC_MIC
% and chain.prior_sd
%
% Plot:
%   left  = QNC Emax
%   right = QNC concentration scale
%   posterior = coloured violin
%   prior = black median + 5/95 percentile interval

%% ========================================================================
% TOGGLES
% ========================================================================

cfg = struct();

cfg.burn_frac = 0.5;
cfg.font_size = 14;
cfg.figure_position = [100 100 1200 550];

cfg.prior_mu = 0;             % log-multiplier prior median
cfg.prior_lo = 5;             % prior lower percentile
cfg.prior_hi = 95;            % prior upper percentile

cfg.violin_width = 0.34;
cfg.face_alpha = 0.28;
cfg.edge_alpha = 0.70;
cfg.violin_line_width = 1.3;

cfg.prior_line_width = 2.0;
cfg.prior_cap_halfwidth = 0.10;
cfg.prior_dot_size = 36;

cfg.show_titles = true;
cfg.show_baseline = true;
cfg.y_lim = [];               % e.g. [-2.5 2.5]

C = lines(max(7, numel(cases)));
for k = 1:numel(cases)
    if ~isfield(cases(k), 'color') || isempty(cases(k).color)
        cases(k).color = C(k,:);
    end
end

%% ========================================================================
% Parameter names
% ========================================================================

param_names = {'log_QNC_Emax', 'log_QNC_MIC'};
param_titles = {'Quinacrine E_{max}', 'Quinacrine concentration scale'};

n_cases = numel(cases);
n_param = numel(param_names);

samples = cell(n_cases, n_param);
prior_sd = NaN(n_cases, n_param);

summary = table();

%% ========================================================================
% Load chains
% ========================================================================

for k = 1:n_cases
    chain = local_load_chain(cases(k).chain_file);

    Theta = chain.theta;
    burn_idx = floor(cfg.burn_frac * size(Theta,1)) + 1;
    Theta = Theta(burn_idx:end,:);

    names = chain.param_names;
    if isstring(names)
        names = cellstr(names);
    end

    for p = 1:n_param
        idx = find(strcmp(names, param_names{p}), 1);
        if isempty(idx)
            error('Parameter %s not found in %s.', ...
                param_names{p}, cases(k).chain_file);
        end

        s = Theta(:,idx);
        samples{k,p} = s;

        if isfield(chain, 'prior_sd') && numel(chain.prior_sd) >= idx
            prior_sd(k,p) = chain.prior_sd(idx);
        else
            error('prior_sd missing or too short in %s.', cases(k).chain_file);
        end

        row = table();
        row.label = string(cases(k).label);
        row.parameter = string(param_names{p});
        row.prior_sd = prior_sd(k,p);
        row.posterior_mean_log = mean(s);
        row.posterior_median_log = median(s);
        row.posterior_q025_log = prctile(s,2.5);
        row.posterior_q25_log = prctile(s,25);
        row.posterior_q75_log = prctile(s,75);
        row.posterior_q975_log = prctile(s,97.5);
        row.posterior_median_multiplier = exp(median(s));

        summary = [summary; row]; %#ok<AGROW>
    end
end

disp(summary)

%% ========================================================================
% Plot
% ========================================================================

fig = figure('Color','w','Position',cfg.figure_position);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for p = 1:n_param
    ax = nexttile;
    hold(ax, 'on')

    for k = 1:n_cases
        s = samples{k,p};

        % posterior violin
        local_violin(ax, ...
            k, s, cases(k).color, ...
            cfg.violin_width, cfg.face_alpha, ...
            cfg.edge_alpha, cfg.violin_line_width);

        % prior interval: median + 5/95 percentiles
        mu = cfg.prior_mu;
        sd = prior_sd(k,p);

        qlo = norminv(cfg.prior_lo/100, mu, sd);
        qhi = norminv(cfg.prior_hi/100, mu, sd);
        qmed = mu;

        % vertical bar
        plot(ax, [k k], [qlo qhi], 'k-', ...
            'LineWidth', cfg.prior_line_width, ...
            'HandleVisibility','off');

        % caps
        w = cfg.prior_cap_halfwidth;
        plot(ax, [k-w k+w], [qlo qlo], 'k-', ...
            'LineWidth', cfg.prior_line_width, ...
            'HandleVisibility','off');
        plot(ax, [k-w k+w], [qhi qhi], 'k-', ...
            'LineWidth', cfg.prior_line_width, ...
            'HandleVisibility','off');

        % median dot
        scatter(ax, k, qmed, cfg.prior_dot_size, ...
            'k', 'filled', ...
            'HandleVisibility','off');
    end

    if cfg.show_baseline
        yline(ax, 0, ':', 'Baseline', ...
            'Color',[0.35 0.35 0.35], ...
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
        title(ax, param_titles{p}, ...
            'FontSize', cfg.font_size + 1, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'tex');
    end

    if ~isempty(cfg.y_lim)
        ylim(ax, cfg.y_lim)
    end
end

% no global title by default

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
    obj = S.(fn{i});
    if isstruct(obj) && isfield(obj, 'theta') && isfield(obj, 'param_names')
        chain = obj;
        return
    end
end

error('No valid chain struct found in %s.', chain_file);

end

%% ========================================================================
function local_violin(ax, x0, s, color, max_width, face_alpha, edge_alpha, line_width)

s = s(:);
s = s(isfinite(s));

if numel(s) < 5
    return
end

lo = prctile(s, 0.5);
hi = prctile(s, 99.5);
pad = 0.15 * max(1e-6, hi-lo);

ygrid = linspace(lo-pad, hi+pad, 250);

try
    f = ksdensity(s, ygrid);
catch
    [f,ygrid] = ksdensity(s);
end

if max(f) > 0
    f = f ./ max(f) * max_width;
end

x_left  = x0 - f;
x_right = x0 + f;

x_patch = [x_left fliplr(x_right)];
y_patch = [ygrid fliplr(ygrid)];

h = patch(ax, x_patch, y_patch, color, ...
    'FaceAlpha', face_alpha, ...
    'EdgeColor', color, ...
    'LineWidth', line_width, ...
    'HandleVisibility','off');

if isprop(h, 'EdgeAlpha')
    h.EdgeAlpha = edge_alpha;
end

end

%{
C = lines(7);

cases = struct([]);

cases(1).chain_file = 'mcmc_sensitivity_outputs/qnc_refractory_51_54_prior_sd0p60_0p60_onechain_mcmc251_vp1001_currentSimulator.mat';
cases(1).label = 'sd = 0.60';
cases(1).color = C(1,:);

cases(2).chain_file = 'mcmc_sensitivity_outputs/qnc_refractory_51_54_prior_sd0p90_0p90_onechain_mcmc252_vp1001_currentSimulator.mat';
cases(2).label = 'sd = 0.90';
cases(2).color = C(2,:);

cases(3).chain_file = 'mcmc_sensitivity_outputs/qnc_refractory_51_54_prior_sd1p25_1p25_onechain_mcmc253_vp1001_currentSimulator.mat';
cases(3).label = 'sd = 1.25';
cases(3).color = C(3,:);

cases(4).chain_file = 'mcmc_sensitivity_outputs/qnc_refractory_51_54_prior_sd1p50_1p50_onechain_mcmc254_vp1001_currentSimulator.mat';
cases(4).label = 'sd = 1.50';
cases(4).color = C(4,:);

summary_prior = plot_qnc_prior_sensitivity_violins(cases);
%}