function summary = plot_project_T_overlay_cases(cases)
%PLOT_PROJECT_T_OVERLAY_CASES
% Overlay posterior predictive T(t) projections from one or more chains.
%
% Each cases(k) must contain:
%   .chain
%   .regimen_name
%   .scenario
%   .absorbExtinction
%   .label
%
% Optional:
%   .line_style

%% ========================================================================
%  USER TOGGLES
%  ========================================================================

font_size = 16;

cfg = struct();

cfg.t_end_days   = 180;
cfg.n_sims       = 800;
cfg.seed         = 123;
cfg.response_day = 30;

cfg.y_scale = 'log';          % 'log' or 'linear'
cfg.show_50_band = true;
cfg.show_95_band = false;
cfg.show_median  = true;

cfg.plot_floor = 1e-16;

cfg.x_lim = [0 cfg.t_end_days];

if strcmpi(cfg.y_scale, 'log')
    cfg.y_lim = [cfg.plot_floor 1e3];
else
    cfg.y_lim = [];
end

cfg.show_clearance_threshold = true;
cfg.show_extinction_threshold = true;
cfg.show_response_day = true;

cfg.extinction_threshold = 1e-9;
cfg.extinction_hold_h    = 72;

%% ========================================================================
%  Validate cases
%  ========================================================================

required = {'chain','regimen_name','scenario','absorbExtinction','label'};

for cc = 1:numel(cases)
    for r = 1:numel(required)
        if ~isfield(cases(cc), required{r})
            error('cases(%d) is missing required field: %s', cc, required{r});
        end
    end

    if ~isfield(cases(cc), 'line_style') || isempty(cases(cc).line_style)
        cases(cc).line_style = '-';
    end
end

rng(cfg.seed);

p0 = giardia_params();

dt_h = 0.5;
if isfield(cases(1).chain, 'scenario') && isfield(cases(1).chain.scenario, 'dt_h')
    dt_h = cases(1).chain.scenario.dt_h;
end

tvec = 0:dt_h:(24*cfg.t_end_days);
time_days = tvec / 24;

%% ========================================================================
%  Simulate each case
%  ========================================================================

case_out = struct([]);

for cc = 1:numel(cases)

    chain = cases(cc).chain;
    regimen_name = cases(cc).regimen_name;
    scenario = cases(cc).scenario;
    absorbExtinction = cases(cc).absorbExtinction;

    % Extract posterior samples
    if isfield(chain, 'theta')
        Theta = chain.theta;
    elseif isfield(chain, 'samples')
        Theta = chain.samples;
    else
        error('Could not find posterior samples in cases(%d).chain.', cc);
    end

    burn = floor(size(Theta,1)/2) + 1;
    Theta = Theta(burn:end,:);

    param_names = chain.param_names;
    SC = chain.scenario;

    % Same VP panel seed for all cases, but panel generation uses each chain's SC.
    % This gives paired virtual-patient draws where the SC structure is comparable.
    Z = local_make_virtual_panel(cfg.n_sims, cfg.seed + 1000);

    Tmat = NaN(cfg.n_sims, numel(tvec));

    success_day = false(cfg.n_sims,1);
    op_clear = false(cfg.n_sims,1);
    rebound = false(cfg.n_sims,1);
    durable = false(cfg.n_sims,1);
    clearance_day = NaN(cfg.n_sims,1);
    rebound_day = NaN(cfg.n_sims,1);
    extinction_day = NaN(cfg.n_sims,1);

    for i = 1:cfg.n_sims

        idx = randi(size(Theta,1));
        theta = Theta(idx,:);

        p = local_make_patient(p0, Z, i, SC);
        p = local_apply_pd_multipliers(p, theta, param_names);

        opts = struct('useMicrobiome', true, 'resistant', false);

        if strcmpi(scenario, 'refractory')
            opts.resistant = true;
        end

        opts.absorbExtinction = absorbExtinction;
        opts.extinction_threshold = cfg.extinction_threshold;
        opts.extinction_hold_h = cfg.extinction_hold_h;

        reg = giardia_regimens(regimen_name, 0);
        out = giardia_simulate(tvec, reg, p, opts);

        T = out.T(:)';
        Tmat(i,:) = T;

        Tday = interp1(time_days, T, cfg.response_day, 'linear', 'extrap');

        cleared_by_day = ~isnan(out.clearance_time_h) && ...
            out.clearance_time_h <= 24*cfg.response_day;

        success_day(i) = cleared_by_day && (Tday < p.clearance_threshold);

        op_clear(i) = ~isnan(out.clearance_time_h);

        if op_clear(i)
            clearance_day(i) = out.clearance_time_h / 24;

            after_idx = find(tvec > out.clearance_time_h);
            if any(T(after_idx) > p.clearance_threshold)
                rebound(i) = true;
                first_rebound_idx = after_idx(find(T(after_idx) > p.clearance_threshold, 1, 'first'));
                rebound_day(i) = tvec(first_rebound_idx) / 24;
            end
        end

        durable(i) = op_clear(i) && ~rebound(i);

        if isfield(out, 'extinction_time_h') && ~isnan(out.extinction_time_h)
            extinction_day(i) = out.extinction_time_h / 24;
        end
    end

    Tplot = Tmat;
    if strcmpi(cfg.y_scale, 'log')
        Tplot = max(Tplot, cfg.plot_floor);
    end

    case_out(cc).Tmat = Tmat;
    case_out(cc).Tplot = Tplot;

    case_out(cc).q025 = prctile(Tplot, 2.5, 1);
    case_out(cc).q25  = prctile(Tplot, 25, 1);
    case_out(cc).q50  = prctile(Tplot, 50, 1);
    case_out(cc).q75  = prctile(Tplot, 75, 1);
    case_out(cc).q975 = prctile(Tplot, 97.5, 1);

    case_out(cc).success_day = success_day;
    case_out(cc).op_clear = op_clear;
    case_out(cc).rebound = rebound;
    case_out(cc).durable = durable;
    case_out(cc).clearance_day = clearance_day;
    case_out(cc).rebound_day = rebound_day;
    case_out(cc).extinction_day = extinction_day;
end

%% ========================================================================
%  Plot
%  ========================================================================

figure('Color','w','Position',[100 100 1150 700]);
hold on

line_width = 2.4;
%{
gray_levels = linspace(0.05, 0.55, max(numel(cases),2));

for cc = 1:numel(cases)

    q025 = case_out(cc).q025;
    q25  = case_out(cc).q25;
    q50  = case_out(cc).q50;
    q75  = case_out(cc).q75;
    q975 = case_out(cc).q975;

    g = gray_levels(cc);
    linecol = [g g g];

    band50 = min(0.85, linecol + 0.35);
    band95 = min(0.92, linecol + 0.55);

    if cfg.show_95_band
        fill([time_days fliplr(time_days)], ...
             [q025 fliplr(q975)], ...
             band95, ...
             'EdgeColor','none', ...
             'FaceAlpha',0.20, ...
             'HandleVisibility','off');
    end

    if cfg.show_50_band
        fill([time_days fliplr(time_days)], ...
             [q25 fliplr(q75)], ...
             band50, ...
             'EdgeColor','none', ...
             'FaceAlpha',0.22, ...
             'HandleVisibility','off');
    end

    if cfg.show_median
        plot(time_days, q50, ...
            'LineStyle', cases(cc).line_style, ...
            'Color', linecol, ...
            'LineWidth', line_width, ...
            'DisplayName', cases(cc).label);
    end
end
%}

if cfg.show_clearance_threshold
    yline(p0.clearance_threshold, ':', 'Operational clearance', ...
        'Color',[0.25 0.25 0.25], ...
        'LineWidth',1.2, ...
        'HandleVisibility','off');
end

if cfg.show_extinction_threshold && strcmpi(cfg.y_scale, 'log')
    yline(cfg.extinction_threshold, ':', 'Extinction threshold', ...
        'Color',[0.45 0.45 0.45], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');
end

if strcmpi(cfg.y_scale, 'log')
    yline(cfg.plot_floor, '-', 'Zero / absorbed shown at floor', ...
        'Color',[0.55 0.55 0.55], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');
end

if cfg.show_response_day
    xline(cfg.response_day, '--', 'Observation day', ...
        'Color',[0.35 0.35 0.35], ...
        'LineWidth',1.1, ...
        'HandleVisibility','off');
end

set(gca, 'YScale', cfg.y_scale, 'fontsize', font_size);

xlabel('Time since treatment start (days)');
ylabel('Trophozoite burden, T(t)');
title('Posterior predictive T(t)', 'Interpreter','none');

grid on
box on
legend('Location','northeast');

if ~isempty(cfg.x_lim)
    xlim(cfg.x_lim);
end

if ~isempty(cfg.y_lim)
    ylim(cfg.y_lim);
end

%% ========================================================================
%  Summary table
%  ========================================================================

summary = table();

Cdefault = lines(max(7, numel(cases)));

for cc = 1:numel(cases)

    q025 = case_out(cc).q025;
    q25  = case_out(cc).q25;
    q50  = case_out(cc).q50;
    q75  = case_out(cc).q75;
    q975 = case_out(cc).q975;

    % Choose colour
    if isfield(cases(cc), 'color') && ~isempty(cases(cc).color)
        thisColor = cases(cc).color;
    else
        thisColor = Cdefault(cc,:);
    end

    % Choose line style
    if isfield(cases(cc), 'line_style') && ~isempty(cases(cc).line_style)
        thisStyle = cases(cc).line_style;
    else
        thisStyle = '-';
    end

    if cfg.show_95_band
        fill([time_days fliplr(time_days)], ...
             [q025 fliplr(q975)], ...
             thisColor, ...
             'EdgeColor','none', ...
             'FaceAlpha',0.10, ...
             'HandleVisibility','off');
    end

    if cfg.show_50_band
        fill([time_days fliplr(time_days)], ...
             [q25 fliplr(q75)], ...
             thisColor, ...
             'EdgeColor','none', ...
             'FaceAlpha',0.18, ...
             'HandleVisibility','off');
    end

    if cfg.show_median
        plot(time_days, q50, ...
            'LineStyle', thisStyle, ...
            'Color', thisColor, ...
            'LineWidth', 3, ... % 2.5, ...
            'DisplayName', cases(cc).label);
    end
end

disp(summary)

end

% -------------------------------------------------------------------------
function Z = local_make_virtual_panel(n_vp, seed)

rng(seed);
Z.rT       = randn(n_vp,1);
Z.K        = randn(n_vp,1);
Z.k_encyst = randn(n_vp,1);
Z.k_c      = randn(n_vp,1);
Z.k_I      = randn(n_vp,1);
Z.k_micro  = randn(n_vp,1);
Z.dys      = randn(n_vp,1);
Z.rM       = randn(n_vp,1);
Z.T0       = randn(n_vp,1);
Z.C0       = randn(n_vp,1);
Z.M0       = randn(n_vp,1);

Z.gut_GE   = randn(n_vp,1);
Z.gut_tr   = randn(n_vp,1);
Z.gut_V    = randn(n_vp,1);
Z.gut_F    = randn(n_vp,1);
Z.alpha    = randn(n_vp,1);

end




%{
S = load('mcmc_outputs/pooled_QNC7D_4chains_currentSimulator.mat');

C = lines(7);

cases = struct([]);

cases(1).chain = S.chain;
cases(1).regimen_name = 'QNC7D';
cases(1).scenario = 'refractory';
cases(1).absorbExtinction = false;
cases(1).label = 'QNC7D, residual persistence';
cases(1).line_style = '-';
cases(1).color = C(1,:);

cases(2).chain = S.chain;
cases(2).regimen_name = 'QNC7D';
cases(2).scenario = 'refractory';
cases(2).absorbExtinction = true;
cases(2).label = 'QNC7D, absorbing extinction';
cases(2).line_style = '-';
cases(2).color = C(7,:);

summary_qnc = plot_project_T_overlay_cases(cases);


S = load('mcmc_outputs/pooled_MTZ14D_ABZ7D_4chains_currentSimulator.mat');

cases = struct([]);

cases(1).chain = S.chain;
cases(1).regimen_name = 'MTZ14D_ABZ7D';
cases(1).scenario = 'refractory';
cases(1).absorbExtinction = false;
cases(1).label = 'MTZ14D+ABZ7D, residual persistence';
cases(1).line_style = '-';
cases(1).color = C(1,:);

cases(2).chain = S.chain;
cases(2).regimen_name = 'MTZ14D_ABZ7D';
cases(2).scenario = 'refractory';
cases(2).absorbExtinction = true;
cases(2).label = 'MTZ14D+ABZ7D, absorbing extinction';
cases(2).line_style = '-';
cases(2).color = C(7,:);

summary_qnc = plot_project_T_overlay_cases(cases);
%}