function summary = plot_project_T_overlay_cases(cases)
%PLOT_PROJECT_T_OVERLAY_CASES
% Overlay posterior predictive T(t) projections from one or more chains.
%
% Each cases(k) must contain:
%   .chain
%   .regimen_name
%   .scenario
%   .absorbExtinction
%   .extinction_threshold
%   .extinction_hold_h
%   .label
%
% Optional:
%   .plot_title
%   .line_style
%   .color
%
% For refractory projections, virtual patients are constructed using the
% same pretreatment framework as giardia_montecarlo:
%
%   candidate VP generation
%       -> 30-day untreated pretreatment burn-in
%       -> retain patients above clearance threshold at treatment start
%       -> cache T(0), C(0), M(0)
%       -> carry infection age into treatment simulation
%
% Posterior samples are taken from the second 50% of chain.theta.
%
% Cases using the same chain/regimen/scenario are paired through the same
% VP seed and posterior-draw seed. Thus residual-persistence and absorbing-
% extinction projections differ only through the extinction assumption.


%% ========================================================================
% USER TOGGLES
% ========================================================================

font_size = 16;

cfg = struct();

cfg.t_end_days   = 180;
cfg.n_sims       = 800;
cfg.seed         = 123;
cfg.response_day = 30;

cfg.pre_treatment_days = 30;

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

%% ========================================================================
% Validate cases
% ========================================================================

required = {'chain','regimen_name','scenario','absorbExtinction', ...
    'extinction_threshold','extinction_hold_h','label'};

for cc = 1:numel(cases)

    for r = 1:numel(required)
        if ~isfield(cases(cc), required{r})
            error('cases(%d) is missing required field: %s', ...
                cc, required{r});
        end
    end

    if ~isfield(cases(cc), 'line_style') || ...
            isempty(cases(cc).line_style)
        cases(cc).line_style = '-';
    end
end


p0 = giardia_params();


%% ========================================================================
% Simulate each case
% ========================================================================

case_out = struct([]);

for cc = 1:numel(cases)

    chain = cases(cc).chain;
    regimen_name = cases(cc).regimen_name;
    scenario = cases(cc).scenario;
    absorbExtinction = cases(cc).absorbExtinction;


    %% --------------------------------------------------------------------
    % Posterior samples: discard first 50% ONCE
    % ---------------------------------------------------------------------

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


    %% --------------------------------------------------------------------
    % Simulation grid
    % ---------------------------------------------------------------------

    dt_h = 0.5;

    if isfield(SC, 'dt_h')
        dt_h = SC.dt_h;
    end

    tvec = 0:dt_h:(24*cfg.t_end_days);
    time_days = tvec / 24;


    %% --------------------------------------------------------------------
    % Construct posterior-predictive VP panel
    %
    % The same seed is deliberately used for each case. Therefore cases
    % with the same biological setup receive the same candidate random
    % numbers and hence the same eligible VP panel.
    % ---------------------------------------------------------------------

    vp_seed = cfg.seed + 1000;

    if strcmpi(scenario, 'refractory')

        P = local_make_eligible_burned_panel( ...
            cfg.n_sims, ...
            vp_seed, ...
            p0, ...
            SC, ...
            regimen_name, ...
            scenario, ...
            cfg.pre_treatment_days);

    elseif strcmpi(scenario, 'baseline')

        Z = local_make_virtual_panel(cfg.n_sims, vp_seed);

        P = repmat(p0, cfg.n_sims, 1);

        for i = 1:cfg.n_sims
            P(i) = local_make_patient(p0, Z, i, SC);
        end

    else

        error('Unknown scenario: %s', scenario);

    end


    %% --------------------------------------------------------------------
    % Posterior projection
    %
    % Reset this RNG identically for every case. Paired cases using the
    % same posterior therefore draw the same theta index for each VP.
    % ---------------------------------------------------------------------

    rng(cfg.seed + 2000);


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

        % Start from already prepared treatment-start VP.
        p = P(i);

        p = local_apply_pd_multipliers( ...
            p, theta, param_names);


        opts = struct();
        opts.useMicrobiome = true;
        opts.resistant = strcmpi(scenario, 'refractory');

        opts.absorbExtinction = absorbExtinction;

        if absorbExtinction
            opts.extinction_threshold = cases(cc).extinction_threshold;
            opts.extinction_hold_h = cases(cc).extinction_hold_h;
        end


        reg = giardia_regimens(regimen_name, 0);

        out = giardia_simulate(tvec, reg, p, opts);


        T = out.T(:)';
        Tmat(i,:) = T;


        %% ---------------------------------------------------------------
        % Clinical response at observation day
        % ---------------------------------------------------------------

        Tday = interp1( ...
            time_days, T, cfg.response_day, ...
            'linear', 'extrap');

        cleared_by_day = ...
            ~isnan(out.clearance_time_h) && ...
            out.clearance_time_h <= 24*cfg.response_day;

        success_day(i) = ...
            cleared_by_day && ...
            (Tday < p.clearance_threshold);


        %% ---------------------------------------------------------------
        % Operational clearance / rebound
        % ---------------------------------------------------------------

        op_clear(i) = ~isnan(out.clearance_time_h);

        if op_clear(i)

            clearance_day(i) = out.clearance_time_h / 24;

            after_idx = find(tvec > out.clearance_time_h);

            if ~isempty(after_idx)

                j = find( ...
                    T(after_idx) >= p.clearance_threshold, ...
                    1, 'first');

                if ~isempty(j)

                    rebound(i) = true;

                    first_rebound_idx = after_idx(j);

                    rebound_day(i) = ...
                        tvec(first_rebound_idx) / 24;

                end
            end
        end


        durable(i) = op_clear(i) && ~rebound(i);


        %% ---------------------------------------------------------------
        % Absorbing extinction
        % ---------------------------------------------------------------

        if isfield(out, 'extinction_time_h') && ...
                ~isnan(out.extinction_time_h)

            extinction_day(i) = ...
                out.extinction_time_h / 24;

        end

    end


    %% --------------------------------------------------------------------
    % Quantiles
    % ---------------------------------------------------------------------

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
% Plot
% ========================================================================

%figure('Color','w','Position',[100 100 1150 700]);
figure('Color','w','Position',[100 100 650 400]);
hold on

Cdefault = lines(max(7, numel(cases)));


for cc = 1:numel(cases)

    q025 = case_out(cc).q025;
    q25  = case_out(cc).q25;
    q50  = case_out(cc).q50;
    q75  = case_out(cc).q75;
    q975 = case_out(cc).q975;


    % Choose colour
    if isfield(cases(cc), 'color') && ...
            ~isempty(cases(cc).color)

        thisColor = cases(cc).color;

    else

        thisColor = Cdefault(cc,:);

    end


    % Choose line style
    if isfield(cases(cc), 'line_style') && ...
            ~isempty(cases(cc).line_style)

        thisStyle = cases(cc).line_style;

    else

        thisStyle = '-';

    end


    % 95% posterior predictive band
    if cfg.show_95_band

        fill( ...
            [time_days fliplr(time_days)], ...
            [q025 fliplr(q975)], ...
            thisColor, ...
            'EdgeColor','none', ...
            'FaceAlpha',0.10, ...
            'HandleVisibility','off');

    end


    % 50% posterior predictive band
    if cfg.show_50_band

        fill( ...
            [time_days fliplr(time_days)], ...
            [q25 fliplr(q75)], ...
            thisColor, ...
            'EdgeColor','none', ...
            'FaceAlpha',0.18, ...
            'HandleVisibility','off');

    end


    % Median
    if cfg.show_median

        plot( ...
            time_days, q50, ...
            'LineStyle', thisStyle, ...
            'Color', thisColor, ...
            'LineWidth',3, ...
            'DisplayName',cases(cc).label);

    end

end


%% ------------------------------------------------------------------------
% Reference lines
% -------------------------------------------------------------------------

ext_idx = find([cases.absorbExtinction], 1, 'first');

if ~isempty(ext_idx)
    extinction_threshold_plot = cases(ext_idx).extinction_threshold;
else
    extinction_threshold_plot = NaN;
end

if cfg.show_clearance_threshold

    yline( ...
        p0.clearance_threshold, ...
        ':', ...
        'Operational clearance', ...
        'Color',[0.25 0.25 0.25], ...
        'LineWidth',1.2, ...
        'HandleVisibility','off');

end


if cfg.show_extinction_threshold && ...
        strcmpi(cfg.y_scale, 'log') && ...
        ~isnan(extinction_threshold_plot)

    yline( ...
        extinction_threshold_plot, ...
        ':', ...
        'Extinction threshold', ...
        'Color',[0.45 0.45 0.45], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');

end


if strcmpi(cfg.y_scale, 'log')

    yline( ...
        cfg.plot_floor, ...
        '-', ...
        'Zero / absorbed shown at floor', ...
        'Color',[0.55 0.55 0.55], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');

end


if cfg.show_response_day

    xline( ...
        cfg.response_day, ...
        '--', ...
        'Observation day', ...
        'Color',[0.35 0.35 0.35], ...
        'LineWidth',1.1, ...
        'HandleVisibility','off', ...
        'LabelVerticalAlignment','middle', ...
        'LabelHorizontalAlignment','left');

end


set(gca, ...
    'YScale',cfg.y_scale, ...
    'FontSize',font_size);

xlabel('Time since treatment start (days)');
ylabel('Trophozoite burden, T(t)');

if isfield(cases(1), 'plot_title') && ~isempty(cases(1).plot_title)

    if ~isempty(ext_idx)
        title_string = sprintf('%s, T_{ext} = 10^{%g}', ...
            cases(1).plot_title, ...
            log10(extinction_threshold_plot));
    else
        title_string = cases(1).plot_title;
    end

    title(title_string, ...
        'FontWeight','bold', ...
        'Interpreter','tex');
end

grid on
box on

legend('Location','northeast');
%legend('Location','east');


if ~isempty(cfg.x_lim)
    xlim(cfg.x_lim);
end

if ~isempty(cfg.y_lim)
    ylim(cfg.y_lim);
end


%% ========================================================================
% Summary table
% ========================================================================

summary = table();

for cc = 1:numel(cases)

    row = table();

    row.label = string(cases(cc).label);
    row.regimen = string(cases(cc).regimen_name);
    row.scenario = string(cases(cc).scenario);
    row.absorbExtinction = cases(cc).absorbExtinction;

    row.n_sims = cfg.n_sims;
    row.response_day = cfg.response_day;

    row.p_response_day = ...
        mean(case_out(cc).success_day);

    row.p_operational_clearance = ...
        mean(case_out(cc).op_clear);

    row.p_rebound = ...
        mean(case_out(cc).rebound);

    row.p_durable = ...
        mean(case_out(cc).durable);

    row.p_absorbed_extinct = ...
        mean(~isnan(case_out(cc).extinction_day));


    row.median_clearance_day = ...
        median( ...
        case_out(cc).clearance_day( ...
        ~isnan(case_out(cc).clearance_day)), ...
        'omitnan');


    row.median_rebound_day = ...
        median( ...
        case_out(cc).rebound_day( ...
        ~isnan(case_out(cc).rebound_day)), ...
        'omitnan');


    row.median_extinction_day = ...
        median( ...
        case_out(cc).extinction_day( ...
        ~isnan(case_out(cc).extinction_day)), ...
        'omitnan');


    summary = [summary; row]; %#ok<AGROW>

end

disp(summary)

end


%% =========================================================================
function P = local_make_eligible_burned_panel( ...
    n_vp, vp_seed, p0, SC, ...
    regimen_name, scenario, pre_treatment_days)
% Generate fixed refractory VPs and perform untreated burn-in once.
%
% This mirrors the current giardia_montecarlo VP construction.

n_candidates = 2 * n_vp;

Zraw = local_make_virtual_panel( ...
    n_candidates, vp_seed);


tvec_burn = ...
    (-24*pre_treatment_days):SC.dt_h:0;


reg = giardia_regimens( ...
    regimen_name, 0);


opts = struct();
opts.useMicrobiome = true;
opts.resistant = strcmpi(scenario, 'refractory');


% Ensure cached patients all share the same structure.
p_template = p0;

p_template.infection_age0_h = ...
    24 * pre_treatment_days;


P_candidate = ...
    repmat(p_template, n_candidates, 1);

eligible = false(n_candidates,1);
T_at_treatment = NaN(n_candidates,1);


for i = 1:n_candidates

    p = local_make_patient( ...
        p0, Zraw, i, SC);

    out = giardia_simulate( ...
        tvec_burn, reg, p, opts);

    T_at_treatment(i) = out.T(end);


    if T_at_treatment(i) >= ...
            p.clearance_threshold

        eligible(i) = true;

        p.T0 = out.T(end);
        p.C0 = out.C(end);
        p.M0 = out.M(end);

        % Preserve infection age when treatment simulation restarts at t=0.
        p.infection_age0_h = ...
            24 * pre_treatment_days;

        P_candidate(i) = p;

    end

end


idx = find(eligible);


if numel(idx) < n_vp

    error( ...
        ['Only %d of %d candidate virtual patients remained infected ' ...
         'after burn-in; need %d. Increase n_candidates.'], ...
        numel(idx), n_candidates, n_vp);

end


idx = idx(1:n_vp);

P = P_candidate(idx);


fprintf('\n');
fprintf('=== Projection refractory VP eligibility: %s ===\n', ...
    regimen_name);

fprintf('Candidate patients:        %d\n', ...
    n_candidates);

fprintf('Eligible after burn-in:    %d\n', ...
    sum(eligible));

fprintf('Retained for projection:   %d\n', ...
    n_vp);

fprintf('Candidate eligibility:     %.3f\n', ...
    mean(eligible));

fprintf('Median eligible T(0):      %.4g\n', ...
    median(T_at_treatment(idx)));

fprintf('Minimum eligible T(0):     %.4g\n', ...
    min(T_at_treatment(idx)));

fprintf('\n');

end


%% =========================================================================
function Z = local_make_virtual_panel(n_vp, seed)
% Common random numbers for VP uncertainty.

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
% Build one VP from fixed random numbers.
% Mirrors the current giardia_montecarlo construction.

p = p0;


%% A. Parasite natural history

if SC.sample_parasite_block

    p.rT = max(1e-6, ...
        p.rT * ...
        exp(SC.vp.parasite.rT * Z.rT(i)));

    p.K = max(1e-6, ...
        p.K * ...
        exp(SC.vp.parasite.K * Z.K(i)));

    p.k_encyst = max(0, ...
        p.k_encyst * ...
        exp(SC.vp.parasite.k_encyst * Z.k_encyst(i)));

    p.k_c_decay = max(0, ...
        p.k_c_decay * ...
        exp(SC.vp.parasite.k_c * Z.k_c(i)));

    p.T0 = max(1e-9, ...
        p.T0 * ...
        exp(SC.vp.parasite.T0 * Z.T0(i)));

    p.C0 = max(0, ...
        p.C0 * ...
        exp(SC.vp.parasite.C0 * Z.C0(i)));

end


%% B. Host / microbiome

if SC.sample_host_block

    p.k_I = max(0, ...
        p.k_I * ...
        exp(SC.vp.host.k_I * Z.k_I(i)));

    p.k_micro = max(0, ...
        p.k_micro * ...
        exp(SC.vp.host.k_micro * Z.k_micro(i)));

    p.dys = max(0, ...
        p.dys * ...
        exp(SC.vp.host.dys * Z.dys(i)));

    p.rM = max(0, ...
        p.rM * ...
        exp(SC.vp.host.rM * Z.rM(i)));

    % Avoid exact M = 0, which is an invariant absorbing boundary.
    M_eps = 1e-6;

    p.M0 = min(1, max(M_eps, ...
        p.M0 + ...
        SC.vp.host.M0_sd * Z.M0(i)));

end


%% C. Shared GI/local exposure multipliers

if SC.sample_gut_block

    drugs = {'MTZ','TDZ','ABZ','NTZ','QNC'};

    for d = 1:numel(drugs)

        nm = drugs{d};


        p.(nm).k_GE = max(1e-4, ...
            p.(nm).k_GE * ...
            exp(SC.vp.gut.k_GE * Z.gut_GE(i)));


        p.(nm).k_tr_duod = max(1e-4, ...
            p.(nm).k_tr_duod * ...
            exp(SC.vp.gut.k_tr_duod * Z.gut_tr(i)));


        p.(nm).V_duod_L = ...
            min(2, max(0.02, ...
            p.(nm).V_duod_L * ...
            exp(SC.vp.gut.V_duod_L * Z.gut_V(i))));


        p.(nm).F_duod = ...
            min(1, max(1e-6, ...
            p.(nm).F_duod * ...
            exp(SC.vp.gut.F_duod * Z.gut_F(i))));

    end


    if isfield(p,'exposure') && ...
            isfield(p.exposure,'alpha_plasma')

        p.exposure.alpha_plasma = ...
            min(0.8, max(0, ...
            p.exposure.alpha_plasma + ...
            SC.vp.gut.alpha_sd * Z.alpha(i)));

    end

end

end


%% =========================================================================
function p = local_apply_pd_multipliers( ...
    p, theta, param_names)
% Apply named posterior log-multiplier PD parameters.

mult = struct();

for j = 1:numel(param_names)
    mult.(param_names{j}) = exp(theta(j));
end


% Nitroimidazole class
if isfield(mult,'log_nitro_Emax')

    p.MTZ.Emax = ...
        p.MTZ.Emax * mult.log_nitro_Emax;

    p.TDZ.Emax = ...
        p.TDZ.Emax * mult.log_nitro_Emax;

end


if isfield(mult,'log_nitro_MIC')

    p.MTZ.MIC = ...
        p.MTZ.MIC * mult.log_nitro_MIC;

    p.TDZ.MIC = ...
        p.TDZ.MIC * mult.log_nitro_MIC;

end


% Albendazole
if isfield(mult,'log_ABZ_Emax')

    p.ABZ.Emax = ...
        p.ABZ.Emax * mult.log_ABZ_Emax;

end


if isfield(mult,'log_ABZ_MIC')

    p.ABZ.MIC = ...
        p.ABZ.MIC * mult.log_ABZ_MIC;

end


% Nitazoxanide
if isfield(mult,'log_NTZ_Emax')

    p.NTZ.Emax = ...
        p.NTZ.Emax * mult.log_NTZ_Emax;

end


if isfield(mult,'log_NTZ_MIC')

    p.NTZ.MIC = ...
        p.NTZ.MIC * mult.log_NTZ_MIC;

end


% Quinacrine
if isfield(mult,'log_QNC_Emax')

    p.QNC.Emax = ...
        p.QNC.Emax * mult.log_QNC_Emax;

end


if isfield(mult,'log_QNC_MIC')

    p.QNC.MIC = ...
        p.QNC.MIC * mult.log_QNC_MIC;

end

end


%{
%% ========================================================================
% HOW TO RUN FROM MATLAB COMMAND WINDOW
% ========================================================================

C = lines(7);


%% ------------------------------------------------------------------------
% QNC7D: residual persistence versus absorbing extinction
% -------------------------------------------------------------------------

clear cases S

S = load( ...
    'mcmc_outputs/pooled_QNC7D_4chains_currentSimulator.mat');


cases = struct([]);


cases(1).chain = S.chain;
cases(1).regimen_name = 'QNC7D';
cases(1).scenario = 'refractory';
cases(1).absorbExtinction = false;
cases(1).extinction_threshold = 1e-12;
cases(1).extinction_hold_h = 72;
cases(1).plot_title = 'QNC7';
cases(1).label = 'QNC7D, residual persistence';
cases(1).line_style = '-';
cases(1).color = C(1,:);


cases(2).chain = S.chain;
cases(2).regimen_name = 'QNC7D';
cases(2).scenario = 'refractory';
cases(2).absorbExtinction = true;
cases(2).extinction_threshold = 1e-12;
cases(2).extinction_hold_h = 72;
cases(2).plot_title = 'QNC7';
cases(2).label = 'QNC7D, absorbing extinction';
cases(2).line_style = '-';
cases(2).color = C(7,:);

cases(1).label = 'Residual persistence';
cases(2).label = 'Absorbing extinction';


summary_qnc = ...
    plot_project_T_overlay_cases(cases);

movegui(gcf,'center');


% Optional:
%
% exportgraphics(gcf, ...
%     'posterior_predictive_T_QNC7D.pdf', ...
%     'ContentType','vector');
%
% exportgraphics(gcf, ...
%     'posterior_predictive_T_QNC7D.png', ...
%     'Resolution',300);
%
% savefig(gcf, ...
%     'posterior_predictive_T_QNC7D.fig');



%% ------------------------------------------------------------------------
% MTZ14D + ABZ7D: residual persistence versus absorbing extinction
% -------------------------------------------------------------------------

clear cases S

S = load( ...
    'mcmc_outputs/pooled_MTZ14D_ABZ7D_4chains_currentSimulator.mat');


cases = struct([]);


cases(1).chain = S.chain;
cases(1).regimen_name = 'MTZ14D_ABZ7D';
cases(1).scenario = 'refractory';
cases(1).absorbExtinction = false;
cases(1).extinction_threshold = 1e-12;
cases(1).extinction_hold_h = 72;
cases(1).plot_title = 'MTZ14+ABZ7';
cases(1).label = 'MTZ14D+ABZ7D, residual persistence';
cases(1).line_style = '-';
cases(1).color = C(1,:);


cases(2).chain = S.chain;
cases(2).regimen_name = 'MTZ14D_ABZ7D';
cases(2).scenario = 'refractory';
cases(2).absorbExtinction = true;
cases(2).extinction_threshold = 1e-12;
cases(2).extinction_hold_h = 72;
cases(2).plot_title = 'MTZ14+ABZ7';
cases(2).label = 'MTZ14D+ABZ7D, absorbing extinction';
cases(2).line_style = '-';
cases(2).color = C(7,:);

cases(1).label = 'Residual persistence';
cases(2).label = 'Absorbing extinction';


summary_mtz_abz = ...
    plot_project_T_overlay_cases(cases);

movegui(gcf,'center');


% Optional:
%
% exportgraphics(gcf, ...
%     'posterior_predictive_T_MTZ14D_ABZ7D.pdf', ...
%     'ContentType','vector');
%
% exportgraphics(gcf, ...
%     'posterior_predictive_T_MTZ14D_ABZ7D.png', ...
%     'Resolution',300);
%
% savefig(gcf, ...
%     'posterior_predictive_T_MTZ14D_ABZ7D.fig');

%}