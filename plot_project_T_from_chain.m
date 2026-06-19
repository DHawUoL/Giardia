function summary = plot_project_T_from_chain(chain, regimen_name, scenario, t_end_days, n_sims, seed, response_day, absorbExtinction)% PLOT_PROJECT_T_FROM_CHAIN
% Posterior predictive trophozoite trajectories from an MCMC chain.
%
% Example:
%   S = load('mcmc_sensitivity_outputs/qnc_refractory_50_54_chain.mat');
%   summary = plot_project_T_from_chain(S.chain, 'QNC7D', 'refractory', 180, 800, 123, 30);
%
% Or if chain is already in workspace:
%   summary = plot_project_T_from_chain(chain_qnc, 'QNC7D', 'refractory', 180, 800, 123, 30);

if nargin < 2 || isempty(regimen_name)
    regimen_name = chain.targets(1).regimen;
end
if nargin < 3 || isempty(scenario)
    scenario = chain.targets(1).scenario;
end
if nargin < 4 || isempty(t_end_days)
    t_end_days = 180;
end
if nargin < 5 || isempty(n_sims)
    n_sims = 800;
end
if nargin < 6 || isempty(seed)
    seed = 123;
end
if nargin < 7 || isempty(response_day)
    if isfield(chain, 'targets') && ~isempty(chain.targets)
        response_day = chain.targets(1).response_day;
    else
        response_day = 30;
    end
end
if nargin < 8 || isempty(absorbExtinction)
    absorbExtinction = false;
end

rng(seed);

% -------------------------------------------------------------------------
% Extract posterior samples
% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
% Simulate posterior predictive trajectories
% -------------------------------------------------------------------------
Tmat = NaN(n_sims, numel(tvec));
success_day = false(n_sims,1);
op_clear = false(n_sims,1);
rebound = false(n_sims,1);
durable = false(n_sims,1);
clearance_day = NaN(n_sims,1);
rebound_day = NaN(n_sims,1);

Z = local_make_virtual_panel(n_sims, seed + 1000);

for i = 1:n_sims

    idx = randi(size(Theta,1));
    theta = Theta(idx,:);

    p = local_make_patient(p0, Z, i, SC);
    p = local_apply_pd_multipliers(p, theta, param_names);

    opts = struct('useMicrobiome', true, 'resistant', false);
    if strcmpi(scenario, 'refractory')
        opts.resistant = true;
    end
    opts.absorbExtinction = absorbExtinction;
    opts.extinction_threshold = 1e-6;
    opts.extinction_hold_h = 72;

    reg = giardia_regimens(regimen_name, 0);
    out = giardia_simulate(tvec, reg, p, opts);

    T = out.T(:)';
    Tmat(i,:) = T;

    Tday = interp1(time_days, T, response_day, 'linear', 'extrap');

    cleared_by_day = ~isnan(out.clearance_time_h) && ...
        out.clearance_time_h <= 24*response_day;

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
end

% -------------------------------------------------------------------------
% Quantiles
% -------------------------------------------------------------------------
q025 = prctile(Tmat, 2.5, 1);
q25  = prctile(Tmat, 25, 1);
q50  = prctile(Tmat, 50, 1);
q75  = prctile(Tmat, 75, 1);
q975 = prctile(Tmat, 97.5, 1);

thr = p0.clearance_threshold;

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
figure('Color','w','Position',[100 100 1050 650]);
hold on

% 95% band
fill([time_days fliplr(time_days)], ...
     [q025 fliplr(q975)], ...
     [0.85 0.85 0.85], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.6, ...
     'DisplayName','95% posterior predictive band');

% 50% band
fill([time_days fliplr(time_days)], ...
     [q25 fliplr(q75)], ...
     [0.65 0.65 0.65], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.6, ...
     'DisplayName','50% posterior predictive band');

plot(time_days, q50, 'k-', 'LineWidth', 2.2, 'DisplayName','Median');

yline(thr, ':', 'Clearance threshold', ...
    'Color',[0.35 0.35 0.35], ...
    'LineWidth',1.2, ...
    'HandleVisibility','off');

xline(response_day, '--', 'Observation day', ...
    'Color',[0.35 0.35 0.35], ...
    'LineWidth',1.2, ...
    'HandleVisibility','off');

set(gca, 'YScale','log');
xlabel('Time since treatment start (days)');
ylabel('Trophozoite burden, T(t)');
title(sprintf('Posterior predictive T(t): %s, %s scenario', regimen_name, scenario), ...
    'Interpreter','none');

grid on
box on
legend('Location','northeast');

% -------------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------------
summary = table();
summary.regimen = string(regimen_name);
summary.scenario = string(scenario);
summary.n_sims = n_sims;
summary.response_day = response_day;
summary.prob_response_day = mean(success_day);
summary.prob_operational_clearance = mean(op_clear);
summary.prob_rebound_given_clearance = mean(rebound(op_clear));
summary.prob_rebound_overall = mean(rebound);
summary.prob_durable_clearance = mean(durable);
summary.median_clearance_day = median(clearance_day(op_clear), 'omitnan');
summary.median_rebound_day = median(rebound_day(rebound), 'omitnan');
summary.final_day = t_end_days;
summary.median_final_T = median(Tmat(:,end), 'omitnan');

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

%{
S = load('mcmc_sensitivity_outputs/qnc_refractory_50_54_chain.mat');
summary_qnc = plot_project_T_from_chain(S.chain, 'QNC7D', 'refractory', 180, 800, 123, 30);

S = load('mcmc_sensitivity_outputs/pooled_alb_nitro_mtz14_abz7_chain.mat');
summary_mtz_abz = plot_project_T_from_chain(S.chain, 'MTZ14D_ABZ7D', 'refractory', 180, 800, 123, 30);
%}