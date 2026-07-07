function summary = plot_project_T_two_panel_treatment_comparison()
%PLOT_PROJECT_T_TWO_PANEL_TREATMENT_COMPARISON
% Two-panel posterior predictive T(t) plot:
%   Panel A: QNC7D
%   Panel B: MTZ14D+ABZ7D
%
% Each panel shows:
%   untreated baseline
%   residual persistence
%   absorbing extinction
%   50% and 95% posterior predictive bands
%
% Requires:
%   giardia_params.m
%   giardia_regimens.m
%   giardia_simulate.m

%% ========================================================================
%  USER TOGGLES
%  ========================================================================

cfg = struct();

% Files
cfg.qnc_file = 'mcmc_outputs/pooled_QNC7D_4chains_currentSimulator.mat';
cfg.mtzabz_file = 'mcmc_outputs/pooled_MTZ14D_ABZ7D_4chains_currentSimulator.mat';

% Simulation
cfg.t_end_days   = 180;
cfg.n_sims       = 800;
cfg.seed         = 123;
cfg.response_day = 30;

% Extinction rule
cfg.extinction_threshold = 1e-9;   % try 1e-9 or 1e-12
cfg.extinction_hold_h    = 72;

% Display
cfg.y_scale = 'log';
cfg.plot_floor = 1e-16;
cfg.y_lim = [cfg.plot_floor 1e3];
cfg.x_lim = [0 cfg.t_end_days];

font_size = 14;
cfg.font_size = font_size;
cfg.title_font_size = font_size + 2;
cfg.line_width = 2.8;

cfg.show_titles = false;       % set true if you want panel titles
cfg.show_50_band = true;
cfg.show_95_band = true;
cfg.show_legend = true;
cfg.show_titles = true;

cfg.figure_position = [100 100 1450 650];

% -------------------------------------------------------------------------
% Colours: edit here
% -------------------------------------------------------------------------
C = lines(7);

col = struct();
col.untreated = [0.15 0.15 0.15];
col.residual  = C(1,:);          % blue
col.extinct   = C(7,:);          % red/brown

% Bands
alpha_95 = 0.06;
alpha_50 = 0.16;

%% ========================================================================
%  Load chains and define panels
%  ========================================================================

Sq = load(cfg.qnc_file, 'chain');
Sm = load(cfg.mtzabz_file, 'chain');

panels = struct([]);

panels(1).chain = Sq.chain;
panels(1).regimen = 'QNC7D';
panels(1).scenario = 'refractory';
panels(1).title = 'QNC7D';

panels(2).chain = Sm.chain;
panels(2).regimen = 'MTZ14D_ABZ7D';
panels(2).scenario = 'refractory';
panels(2).title = 'MTZ14D+ABZ7D';

p0 = giardia_params();

dt_h = 0.5;
if isfield(Sq.chain, 'scenario') && isfield(Sq.chain.scenario, 'dt_h')
    dt_h = Sq.chain.scenario.dt_h;
end

tvec = 0:dt_h:(24*cfg.t_end_days);
time_days = tvec / 24;

summary = table();

%% ========================================================================
%  Figure
%  ========================================================================

fig = figure('Color','w','Position',cfg.figure_position);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for pp = 1:numel(panels)

    ax = nexttile;
    hold(ax, 'on')

    chain = panels(pp).chain;

    % ---------------------------------------------------------------------
    % Untreated baseline
    % ---------------------------------------------------------------------
    untreated = local_project_distribution( ...
        [], ...
        chain, ...
        'UNTREATED', ...
        panels(pp).scenario, ...
        false, ...
        cfg, ...
        tvec, ...
        cfg.seed + 5000 + pp);

    cfg_no_band = cfg;
        cfg_no_band.show_50_band = false;
        cfg_no_band.show_95_band = false;
        
        local_plot_band_and_median( ...
            ax, time_days, untreated, ...
            col.untreated, ...
            '--', ...
            'Untreated baseline', ...
            cfg_no_band, alpha_95, alpha_50, true);

    % ---------------------------------------------------------------------
    % Treated: residual persistence
    % ---------------------------------------------------------------------
    residual = local_project_distribution( ...
        chain, ...
        chain, ...
        panels(pp).regimen, ...
        panels(pp).scenario, ...
        false, ...
        cfg, ...
        tvec, ...
        cfg.seed + 1000 + pp);

    local_plot_band_and_median( ...
        ax, time_days, residual, ...
        col.residual, ...
        '-', ...
        'Residual persistence', ...
        cfg, alpha_95, alpha_50, true);

    % ---------------------------------------------------------------------
    % Treated: absorbing extinction
    % ---------------------------------------------------------------------
    extinct = local_project_distribution( ...
        chain, ...
        chain, ...
        panels(pp).regimen, ...
        panels(pp).scenario, ...
        true, ...
        cfg, ...
        tvec, ...
        cfg.seed + 2000 + pp);

    local_plot_band_and_median( ...
        ax, time_days, extinct, ...
        col.extinct, ...
        '-', ...
        sprintf('Absorbing extinction, T_{ext}=%.0e', cfg.extinction_threshold), ...
        cfg, alpha_95, alpha_50, true);

    % ---------------------------------------------------------------------
    % Reference lines
    % ---------------------------------------------------------------------
    yline(ax, p0.clearance_threshold, ':', 'Operational clearance', ...
        'Color',[0.25 0.25 0.25], ...
        'LineWidth',1.2, ...
        'HandleVisibility','off');

    yline(ax, cfg.extinction_threshold, ':', 'Extinction threshold', ...
        'Color',[0.45 0.45 0.45], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');

    yline(ax, cfg.plot_floor, '-', ...
        'Color',[0.55 0.55 0.55], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');

    xline(ax, cfg.response_day, '--', 'Observation day', ...
        'Color',[0.35 0.35 0.35], ...
        'LineWidth',1.1, ...
        'HandleVisibility','off');

    % ---------------------------------------------------------------------
    % Axes
    % ---------------------------------------------------------------------
    set(ax, 'YScale', cfg.y_scale, 'FontSize', cfg.font_size);
    xlim(ax, cfg.x_lim);
    ylim(ax, cfg.y_lim);

    xlabel(ax, 'Time since treatment start (days)');
    ylabel(ax, 'Trophozoite burden, T(t)');

    grid(ax, 'on')
    box(ax, 'on')

    if cfg.show_titles
        title(ax, panels(pp).title, ...
        'FontSize', cfg.title_font_size, ...
        'FontWeight','bold', ...
        'Interpreter','none');
    end

    if cfg.show_legend
        legend(ax, 'Location','southeast', 'FontSize', cfg.font_size - 1);
    end

    % ---------------------------------------------------------------------
    % Summary rows
    % ---------------------------------------------------------------------
    summary = [summary; local_summary_row(panels(pp).title, 'Untreated baseline', untreated, cfg, p0)]; %#ok<AGROW>
    summary = [summary; local_summary_row(panels(pp).title, 'Residual persistence', residual, cfg, p0)]; %#ok<AGROW>
    summary = [summary; local_summary_row(panels(pp).title, 'Absorbing extinction', extinct, cfg, p0)]; %#ok<AGROW>
end

%sgtitle(tl, 'Posterior predictive T(t)', ...
%    'FontWeight','bold', ...
%    'FontSize', cfg.title_font_size);

disp(summary)

end

%% =========================================================================
function outdist = local_project_distribution(chain_for_theta, chain_for_SC, regimen_name, scenario, absorbExtinction, cfg, tvec, seed)

rng(seed);

p0 = giardia_params();

time_days = tvec / 24;
n_t = numel(tvec);

SC = chain_for_SC.scenario;

% Posterior samples
if isempty(chain_for_theta)
    Theta = [];
    param_names = {};
else
    if isfield(chain_for_theta, 'theta')
        Theta = chain_for_theta.theta;
    elseif isfield(chain_for_theta, 'samples')
        Theta = chain_for_theta.samples;
    else
        error('Could not find posterior samples in chain.');
    end

    burn = floor(size(Theta,1)/2) + 1;
    Theta = Theta(burn:end,:);
    param_names = chain_for_theta.param_names;
end

Z = local_make_virtual_panel(cfg.n_sims, seed + 1000);

Tmat = NaN(cfg.n_sims, n_t);

clearance_time_h = NaN(cfg.n_sims,1);
rebound_time_h = NaN(cfg.n_sims,1);
extinction_time_h = NaN(cfg.n_sims,1);

for i = 1:cfg.n_sims

    p = local_make_patient(p0, Z, i, SC);

    if ~isempty(chain_for_theta)
        idx = randi(size(Theta,1));
        theta = Theta(idx,:);
        p = local_apply_pd_multipliers(p, theta, param_names);
    end

    opts = struct('useMicrobiome', true, 'resistant', false);

    if strcmpi(scenario, 'refractory')
        opts.resistant = true;
    elseif strcmpi(scenario, 'baseline')
        opts.resistant = false;
    else
        error('Unknown scenario: %s', scenario);
    end

    opts.absorbExtinction = absorbExtinction;
    opts.extinction_threshold = cfg.extinction_threshold;
    opts.extinction_hold_h = cfg.extinction_hold_h;

    if strcmpi(regimen_name, 'UNTREATED')
        reg = struct();
        reg.name = 'UNTREATED';
        reg.doses = struct();
    else
        reg = giardia_regimens(regimen_name, 0);
    end

    out = giardia_simulate(tvec, reg, p, opts);

    T = out.T(:)';
    Tmat(i,:) = T;

    if isfield(out, 'clearance_time_h') && ~isnan(out.clearance_time_h)
        clearance_time_h(i) = out.clearance_time_h;

        after_idx = find(tvec > out.clearance_time_h);
        rb = after_idx(find(T(after_idx) >= p.clearance_threshold, 1, 'first'));

        if ~isempty(rb)
            rebound_time_h(i) = tvec(rb);
        end
    end

    if isfield(out, 'extinction_time_h') && ~isnan(out.extinction_time_h)
        extinction_time_h(i) = out.extinction_time_h;
    end
end

Tplot = Tmat;
if strcmpi(cfg.y_scale, 'log')
    Tplot = max(Tplot, cfg.plot_floor);
end

outdist = struct();
outdist.time_days = time_days;
outdist.Tmat = Tmat;
outdist.Tplot = Tplot;

outdist.q025 = prctile(Tplot, 2.5, 1);
outdist.q25  = prctile(Tplot, 25, 1);
outdist.q50  = prctile(Tplot, 50, 1);
outdist.q75  = prctile(Tplot, 75, 1);
outdist.q975 = prctile(Tplot, 97.5, 1);

outdist.clearance_time_h = clearance_time_h;
outdist.rebound_time_h = rebound_time_h;
outdist.extinction_time_h = extinction_time_h;

end

%% =========================================================================
function local_plot_band_and_median(ax, time_days, D, color, line_style, label, cfg, alpha_95, alpha_50, show_in_legend)

if cfg.show_95_band
    fill(ax, [time_days fliplr(time_days)], ...
        [D.q025 fliplr(D.q975)], ...
        color, ...
        'EdgeColor','none', ...
        'FaceAlpha',alpha_95, ...
        'HandleVisibility','off');
end

if cfg.show_50_band
    fill(ax, [time_days fliplr(time_days)], ...
        [D.q25 fliplr(D.q75)], ...
        color, ...
        'EdgeColor','none', ...
        'FaceAlpha',alpha_50, ...
        'HandleVisibility','off');
end

if show_in_legend
    plot(ax, time_days, D.q50, ...
        'LineStyle', line_style, ...
        'Color', color, ...
        'LineWidth', cfg.line_width, ...
        'DisplayName', label);
else
    plot(ax, time_days, D.q50, ...
        'LineStyle', line_style, ...
        'Color', color, ...
        'LineWidth', cfg.line_width, ...
        'HandleVisibility','off');
end

end

%% =========================================================================
function row = local_summary_row(panel_label, assumption_label, D, cfg, p0)

final_h = 24*cfg.t_end_days;

cleared_by_final = ~isnan(D.clearance_time_h) & D.clearance_time_h <= final_h;
rebound_by_final = ~isnan(D.rebound_time_h) & D.rebound_time_h <= final_h;
durable_by_final = cleared_by_final & ~rebound_by_final;

row = table();
row.panel = string(panel_label);
row.assumption = string(assumption_label);
row.n_sims = cfg.n_sims;
row.final_day = cfg.t_end_days;
row.p_operational_clearance_by_final = mean(cleared_by_final);
row.p_rebound_by_final = mean(rebound_by_final);
row.p_durable_clearance_by_final = mean(durable_by_final);
row.p_absorbed_extinct_by_final = mean(~isnan(D.extinction_time_h) & D.extinction_time_h <= final_h);
row.median_final_T = median(D.Tmat(:,end), 'omitnan');
row.operational_clearance_threshold = p0.clearance_threshold;
row.extinction_threshold = cfg.extinction_threshold;

end

%% =========================================================================
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

%% =========================================================================
function p = local_make_patient(p0, Z, i, SC)

p = p0;

if SC.sample_parasite_block
    p.rT        = max(1e-6, p.rT        * exp(SC.vp.parasite.rT       * Z.rT(i)));
    p.K         = max(1e-6, p.K         * exp(SC.vp.parasite.K        * Z.K(i)));
    p.k_encyst  = max(0,    p.k_encyst  * exp(SC.vp.parasite.k_encyst * Z.k_encyst(i)));
    p.k_c_decay = max(0,    p.k_c_decay * exp(SC.vp.parasite.k_c      * Z.k_c(i)));
    p.T0        = max(1e-9, p.T0        * exp(SC.vp.parasite.T0       * Z.T0(i)));
    p.C0        = max(0,    p.C0        * exp(SC.vp.parasite.C0       * Z.C0(i)));
end

if SC.sample_host_block
    p.k_I     = max(0, p.k_I     * exp(SC.vp.host.k_I     * Z.k_I(i)));
    p.k_micro = max(0, p.k_micro * exp(SC.vp.host.k_micro * Z.k_micro(i)));
    p.dys     = max(0, p.dys     * exp(SC.vp.host.dys     * Z.dys(i)));
    p.rM      = max(0, p.rM      * exp(SC.vp.host.rM      * Z.rM(i)));
    p.M0      = min(1, max(0, p.M0 + SC.vp.host.M0_sd * Z.M0(i)));
end

if SC.sample_gut_block
    drugs = {'MTZ','TDZ','ABZ','NTZ','QNC'};
    for d = 1:numel(drugs)
        nm = drugs{d};
        p.(nm).k_GE      = max(1e-4, p.(nm).k_GE      * exp(SC.vp.gut.k_GE      * Z.gut_GE(i)));
        p.(nm).k_tr_duod = max(1e-4, p.(nm).k_tr_duod * exp(SC.vp.gut.k_tr_duod * Z.gut_tr(i)));
        p.(nm).V_duod_L  = min(2, max(0.02, p.(nm).V_duod_L * exp(SC.vp.gut.V_duod_L * Z.gut_V(i))));
        p.(nm).F_duod    = min(1, max(1e-6, p.(nm).F_duod   * exp(SC.vp.gut.F_duod   * Z.gut_F(i))));
    end

    if isfield(p, 'exposure') && isfield(p.exposure, 'alpha_plasma')
        p.exposure.alpha_plasma = min(0.8, max(0, p.exposure.alpha_plasma + SC.vp.gut.alpha_sd * Z.alpha(i)));
    end
end

end

%% =========================================================================
function p = local_apply_pd_multipliers(p, theta, param_names)

mult = struct();
for j = 1:numel(param_names)
    mult.(param_names{j}) = exp(theta(j));
end

if isfield(mult,'log_nitro_Emax')
    p.MTZ.Emax = p.MTZ.Emax * mult.log_nitro_Emax;
    p.TDZ.Emax = p.TDZ.Emax * mult.log_nitro_Emax;
end

if isfield(mult,'log_nitro_MIC')
    p.MTZ.MIC = p.MTZ.MIC * mult.log_nitro_MIC;
    p.TDZ.MIC = p.TDZ.MIC * mult.log_nitro_MIC;
end

if isfield(mult,'log_ABZ_Emax')
    p.ABZ.Emax = p.ABZ.Emax * mult.log_ABZ_Emax;
end

if isfield(mult,'log_ABZ_MIC')
    p.ABZ.MIC = p.ABZ.MIC * mult.log_ABZ_MIC;
end

if isfield(mult,'log_NTZ_Emax')
    p.NTZ.Emax = p.NTZ.Emax * mult.log_NTZ_Emax;
end

if isfield(mult,'log_NTZ_MIC')
    p.NTZ.MIC = p.NTZ.MIC * mult.log_NTZ_MIC;
end

if isfield(mult,'log_QNC_Emax')
    p.QNC.Emax = p.QNC.Emax * mult.log_QNC_Emax;
end

if isfield(mult,'log_QNC_MIC')
    p.QNC.MIC = p.QNC.MIC * mult.log_QNC_MIC;
end

end