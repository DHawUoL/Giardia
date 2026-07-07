function summary = plot_outcome_probability_curves_cases(cases)
%PLOT_OUTCOME_PROBABILITY_CURVES_CASES
% Posterior predictive outcome probabilities for multiple projection cases.
%
% Each cases(k) needs:
%   .chain_file
%   .regimen
%   .scenario
%   .label
%
% Optional per case:
%   .absorbExtinction
%   .extinction_threshold
%   .extinction_hold_h
%   .color
%   .line_style
%
% This is intended for comparing:
%   - QNC7D vs MTZ14D_ABZ7D
%   - residual persistence vs absorbing extinction
%   - extinction thresholds such as 1e-9 and 1e-12

%% ========================================================================
%  USER TOGGLES
%  ========================================================================

cfg = struct();

% Simulation
cfg.t_end_days   = 180;
cfg.n_sims       = 800;
cfg.seed         = 123;
cfg.response_day = 30;

% Default extinction rule if not specified in cases(k)
cfg.default_absorbExtinction = false;
cfg.default_extinction_threshold = NaN;
cfg.default_extinction_hold_h = 72;

% Display
cfg.figure_position = [100 100 1350 850];
cfg.line_width = 4;%2.4;
cfg.font_size = 16;
cfg.show_legend = true;
cfg.legend_num_columns = 1;

% -------------------------------------------------------------------------
% COLOURS: change these here
% -------------------------------------------------------------------------
C = lines(7);

col = struct();
col.QNC       = C(1,:);          % blue
col.MTZABZ    = C(7,:);          % dark red / brown
col.grey      = [0.45 0.45 0.45];
col.black     = [0 0 0];

% Optional threshold-specific variants if desired
col.QNC_resid    = col.QNC;
col.QNC_ext9     = col.QNC;
col.QNC_ext12    = col.QNC;

col.MTZABZ_resid = col.MTZABZ;
col.MTZABZ_ext9  = col.MTZABZ;
col.MTZABZ_ext12 = col.MTZABZ;

% Line styles: change here
ls = struct();
ls.residual = '-';
ls.ext9     = '--';
ls.ext12    = ':';

% Reference lines
cfg.show_response_day = true;
cfg.response_day_color = [0.45 0.45 0.45];

%% ========================================================================
%  Fill missing case fields
%  ========================================================================

for k = 1:numel(cases)

    if ~isfield(cases(k), 'absorbExtinction') || isempty(cases(k).absorbExtinction)
        cases(k).absorbExtinction = cfg.default_absorbExtinction;
    end

    if ~isfield(cases(k), 'extinction_threshold') || isempty(cases(k).extinction_threshold)
        cases(k).extinction_threshold = cfg.default_extinction_threshold;
    end

    if ~isfield(cases(k), 'extinction_hold_h') || isempty(cases(k).extinction_hold_h)
        cases(k).extinction_hold_h = cfg.default_extinction_hold_h;
    end

    if ~isfield(cases(k), 'line_style') || isempty(cases(k).line_style)
        cases(k).line_style = '-';
    end

    if ~isfield(cases(k), 'color') || isempty(cases(k).color)
        cases(k).color = C(mod(k-1, size(C,1)) + 1, :);
    end
end

%% ========================================================================
%  Project each case
%  ========================================================================

n_cases = numel(cases);
all_curves = cell(n_cases,1);
summary = table();

for k = 1:n_cases

    S = load(cases(k).chain_file, 'chain');
    chain = S.chain;

    [curves, row] = local_project_one_chain( ...
        chain, ...
        cases(k).regimen, ...
        cases(k).scenario, ...
        cases(k).label, ...
        cfg.t_end_days, ...
        cfg.n_sims, ...
        cfg.seed + 100*k, ...
        cfg.response_day, ...
        cases(k).absorbExtinction, ...
        cases(k).extinction_threshold, ...
        cases(k).extinction_hold_h);

    all_curves{k} = curves;
    summary = [summary; row]; %#ok<AGROW>
end

disp(summary)
writetable(summary, 'posterior_predictive_outcome_summary_cases.csv');

%% ========================================================================
%  Plot four outcome curves
%  ========================================================================

time_days = all_curves{1}.time_days;

fig = figure('Color','w','Position',cfg.figure_position);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% -------------------------------------------------------------------------
% 1. Below threshold at time t
% -------------------------------------------------------------------------
nexttile; hold on
for k = 1:n_cases
    plot(time_days, all_curves{k}.p_below_now, ...
        'LineWidth', cfg.line_width, ...
        'Color', cases(k).color, ...
        'LineStyle', cases(k).line_style, ...
        'DisplayName', cases(k).label);
end
local_obs_line(cfg)
xlim([0,time_days(end)])
ylim([0 1])
xlabel('Time since treatment start (days)')
ylabel('Probability')
title('Below operational clearance threshold at time t')
grid on; box on

% -------------------------------------------------------------------------
% 2. Operational clearance by time t
% -------------------------------------------------------------------------
nexttile; hold on
for k = 1:n_cases
    plot(time_days, all_curves{k}.p_operational_clearance_by, ...
        'LineWidth', cfg.line_width, ...
        'Color', cases(k).color, ...
        'LineStyle', cases(k).line_style, ...
        'DisplayName', cases(k).label);
end
local_obs_line(cfg)
xlim([0,time_days(end)])
ylim([0 1])
xlabel('Time since treatment start (days)')
ylabel('Probability')
title('Operational clearance by time t')
grid on; box on

% -------------------------------------------------------------------------
% 3. Rebound by time t
% -------------------------------------------------------------------------
nexttile; hold on
for k = 1:n_cases
    plot(time_days, all_curves{k}.p_rebound_by, ...
        'LineWidth', cfg.line_width, ...
        'Color', cases(k).color, ...
        'LineStyle', cases(k).line_style, ...
        'DisplayName', cases(k).label);
end
local_obs_line(cfg)
xlim([0,time_days(end)])
ylim([0 1])
xlabel('Time since treatment start (days)')
ylabel('Probability')
title('Rebound after operational clearance by time t')
grid on; box on

% -------------------------------------------------------------------------
% 4. Durable clearance by time t
% -------------------------------------------------------------------------
nexttile; hold on
for k = 1:n_cases
    plot(time_days, all_curves{k}.p_durable_by, ...
        'LineWidth', cfg.line_width, ...
        'Color', cases(k).color, ...
        'LineStyle', cases(k).line_style, ...
        'DisplayName', cases(k).label);
end
local_obs_line(cfg)
xlim([0,time_days(end)])
ylim([0 1])
xlabel('Time since treatment start (days)')
ylabel('Probability')
title('Durable clearance by time t')
grid on; box on

ax = findall(fig, 'Type', 'axes');
set(ax, 'FontSize', cfg.font_size);

if cfg.show_legend
    lgd = legend('Location','eastoutside', ...
        'NumColumns', cfg.legend_num_columns);
    lgd.Layout.Tile = 'east';
    lgd.FontSize = cfg.font_size;
end

%sgtitle(tl, 'Posterior predictive outcome probabilities', ...
    %'FontWeight','bold', ...
    %'FontSize', cfg.font_size + 2);

end

%% =========================================================================
function local_obs_line(cfg)

if cfg.show_response_day
    xline(cfg.response_day, '--', 'Observation day', ...
    'Color', cfg.response_day_color, ...
    'HandleVisibility','off', ...
    'LabelVerticalAlignment','middle', ...
    'LabelHorizontalAlignment','left');
end

end

%% =========================================================================
function [curves, summary_row] = local_project_one_chain(chain, regimen_name, scenario, label, ...
    t_end_days, n_sims, seed, response_day, absorbExtinction, extinction_threshold, extinction_hold_h)

rng(seed);

% Extract posterior samples
if isfield(chain, 'theta')
    Theta = chain.theta;
elseif isfield(chain, 'samples')
    Theta = chain.samples;
else
    error('Could not find posterior samples in chain.');
end

burn = floor(size(Theta,1)/2) + 1;
Theta = Theta(burn:end,:);

param_names = chain.param_names;
SC = chain.scenario;

p0 = giardia_params();

dt_h = 0.5;
if isfield(SC, 'dt_h')
    dt_h = SC.dt_h;
end

tvec = 0:dt_h:(24*t_end_days);
time_days = tvec / 24;
n_t = numel(tvec);

% Fixed/generated virtual patient panel for this projection
Z = local_make_virtual_panel(n_sims, seed + 1000);

below_now = false(n_sims, n_t);
clearance_time_h = NaN(n_sims,1);
rebound_time_h   = NaN(n_sims,1);
extinction_time_h = NaN(n_sims,1);
response_day_success = false(n_sims,1);

for i = 1:n_sims

    theta = Theta(randi(size(Theta,1)),:);

    p = local_make_patient(p0, Z, i, SC);
    p = local_apply_pd_multipliers(p, theta, param_names);

    opts = struct('useMicrobiome', true, 'resistant', false);

    if strcmpi(scenario, 'refractory')
        opts.resistant = true;
    elseif strcmpi(scenario, 'baseline')
        opts.resistant = false;
    else
        error('Unknown scenario: %s', scenario);
    end

    opts.absorbExtinction = absorbExtinction;

    if absorbExtinction
        opts.extinction_threshold = extinction_threshold;
        opts.extinction_hold_h = extinction_hold_h;
    end

    reg = giardia_regimens(regimen_name, 0);
    out = giardia_simulate(tvec, reg, p, opts);

    T = out.T(:)';
    thr = p.clearance_threshold;

    below_now(i,:) = T < thr;

    if isfield(out, 'clearance_time_h') && ~isnan(out.clearance_time_h)
        clearance_time_h(i) = out.clearance_time_h;

        after_idx = find(tvec > out.clearance_time_h);
        rb = after_idx(find(T(after_idx) >= thr, 1, 'first'));

        if ~isempty(rb)
            rebound_time_h(i) = tvec(rb);
        end
    end

    if isfield(out, 'extinction_time_h') && ~isnan(out.extinction_time_h)
        extinction_time_h(i) = out.extinction_time_h;
    end

    Tday = interp1(time_days, T, response_day, 'linear', 'extrap');
    cleared_by_day = ~isnan(clearance_time_h(i)) && clearance_time_h(i) <= 24*response_day;
    response_day_success(i) = cleared_by_day && (Tday < thr);
end

% Time-dependent outcome probabilities
p_below_now = mean(below_now, 1);

p_operational_clearance_by = NaN(1,n_t);
p_rebound_by = NaN(1,n_t);
p_durable_by = NaN(1,n_t);

for tt = 1:n_t

    th = tvec(tt);

    cleared_by_t = ~isnan(clearance_time_h) & clearance_time_h <= th;
    rebound_by_t = ~isnan(rebound_time_h) & rebound_time_h <= th;

    durable_by_t = cleared_by_t & ~rebound_by_t;

    p_operational_clearance_by(tt) = mean(cleared_by_t);
    p_rebound_by(tt) = mean(rebound_by_t);
    p_durable_by(tt) = mean(durable_by_t);
end

curves = struct();
curves.label = label;
curves.regimen = regimen_name;
curves.scenario = scenario;
curves.absorbExtinction = absorbExtinction;
curves.extinction_threshold = extinction_threshold;
curves.extinction_hold_h = extinction_hold_h;
curves.time_days = time_days;
curves.p_below_now = p_below_now;
curves.p_operational_clearance_by = p_operational_clearance_by;
curves.p_rebound_by = p_rebound_by;
curves.p_durable_by = p_durable_by;

% Summary at observation day and final time
obs_h = 24*response_day;
final_h = 24*t_end_days;

cleared_by_obs = ~isnan(clearance_time_h) & clearance_time_h <= obs_h;
rebound_by_obs = ~isnan(rebound_time_h) & rebound_time_h <= obs_h;
durable_by_obs = cleared_by_obs & ~rebound_by_obs;

cleared_by_final = ~isnan(clearance_time_h) & clearance_time_h <= final_h;
rebound_by_final = ~isnan(rebound_time_h) & rebound_time_h <= final_h;
durable_by_final = cleared_by_final & ~rebound_by_final;

summary_row = table();
summary_row.label = string(label);
summary_row.regimen = string(regimen_name);
summary_row.scenario = string(scenario);
summary_row.absorbExtinction = absorbExtinction;
summary_row.extinction_threshold = extinction_threshold;
summary_row.extinction_hold_h = extinction_hold_h;
summary_row.n_sims = n_sims;
summary_row.response_day = response_day;

summary_row.p_response_day = mean(response_day_success);
summary_row.p_operational_clearance_by_response_day = mean(cleared_by_obs);
summary_row.p_rebound_by_response_day = mean(rebound_by_obs);
summary_row.p_durable_by_response_day = mean(durable_by_obs);

summary_row.final_day = t_end_days;
summary_row.p_operational_clearance_by_final = mean(cleared_by_final);
summary_row.p_rebound_by_final = mean(rebound_by_final);
summary_row.p_durable_by_final = mean(durable_by_final);
summary_row.p_absorbed_extinct_by_final = mean(~isnan(extinction_time_h) & extinction_time_h <= final_h);

summary_row.median_clearance_day = median(clearance_time_h(~isnan(clearance_time_h))/24, 'omitnan');
summary_row.median_rebound_day = median(rebound_time_h(~isnan(rebound_time_h))/24, 'omitnan');
summary_row.median_extinction_day = median(extinction_time_h(~isnan(extinction_time_h))/24, 'omitnan');

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

% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
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