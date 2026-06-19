function chain = giardia_montecarlo(n_iter, n_vp, t_end_days, mcmc_seed, scenario_name, vp_seed, theta_init)% GIARDIA_MONTECARLO
% Random-walk Metropolis calibration of Giardia pharmacodynamic parameters.
%
% Core idea:
%   - PK/exposure and baseline biology are fixed or externally specified.
%   - host/parasite/gut uncertainty is represented by a fixed virtual-patient panel.
%   - PD uncertainty is explored through MCMC log-multipliers.
%
% Usage:
%   chain = giardia_montecarlo();
%   chain = giardia_montecarlo(3000, 150, 35, 42, 'pooled_alb_nitro_mtz14_abz7', 1001);
%
% Arguments:
%   n_iter        number of MCMC iterations
%   n_vp          number of virtual patients
%   t_end_days    simulation horizon used for likelihood target
%   mcmc_seed     seed for MCMC proposals / accept-reject randomness
%   scenario_name calibration scenario name
%   vp_seed       seed for fixed virtual-patient panel
%
% Requirements:
%   giardia_params.m
%   giardia_regimens.m
%   giardia_simulate.m

if nargin < 1 || isempty(n_iter),        n_iter = 2000; end
if nargin < 2 || isempty(n_vp),          n_vp = 120; end
if nargin < 3 || isempty(t_end_days),    t_end_days = 35; end
if nargin < 4 || isempty(mcmc_seed),     mcmc_seed = 42; end
if nargin < 5 || isempty(scenario_name), scenario_name = 'pooled_alb_nitro_mtz14_abz7'; end
if nargin < 6 || isempty(vp_seed),       vp_seed = 1001; end
if nargin < 7 || isempty(theta_init); theta_init = []; end

%% ========================================================================
%  SCENARIO / PARAMETER BLOCK
%  ========================================================================

SC = struct();

switch lower(scenario_name)

    case 'mtz_abz_refractory'
        SC.name = 'MTZ_ABZ_refractory';
        SC.param_names = {'log_nitro_Emax','log_nitro_MIC', ...
                          'log_ABZ_Emax','log_ABZ_MIC'};
        SC.prior_sd = [0.75 0.75 0.50 0.50];
        SC.prop_sd  = [0.24 0.24 0.18 0.18];

        SC.targets = struct([]);
        SC.targets(1).name         = 'Meltzer2014_pooled_Alb_nitro';
        SC.targets(1).regimen      = 'MTZ14D_ABZ7D';
        SC.targets(1).scenario     = 'refractory';
        SC.targets(1).x            = 42;
        SC.targets(1).n            = 53;
        SC.targets(1).response_day = 30;
        SC.targets(1).weight       = 1.0;

    case 'nitro_baseline'
        SC.name = 'nitro_baseline';
        SC.param_names = {'log_nitro_Emax','log_nitro_MIC'};
        SC.prior_sd = [0.75 0.75];
        SC.prop_sd  = [0.24 0.24];

        SC.targets = struct([]);
        SC.targets(1).name         = 'Karabay_MTZ5D_day15';
        SC.targets(1).regimen      = 'MTZ5D';
        SC.targets(1).scenario     = 'baseline';
        SC.targets(1).x            = 29;
        SC.targets(1).n            = 29;
        SC.targets(1).response_day = 15;
        SC.targets(1).weight       = 1.0;

    case 'abz_baseline'
        SC.name = 'ABZ_baseline';
        SC.param_names = {'log_ABZ_Emax','log_ABZ_MIC'};
        SC.prior_sd = [0.50 0.50];
        SC.prop_sd  = [0.18 0.18];

        SC.targets = struct([]);
        SC.targets(1).name         = 'Karabay_ABZ5D_day15';
        SC.targets(1).regimen      = 'ABZ5D';
        SC.targets(1).scenario     = 'baseline';
        SC.targets(1).x            = 27;
        SC.targets(1).n            = 28;
        SC.targets(1).response_day = 15;
        SC.targets(1).weight       = 1.0;

    case 'qnc_refractory'
        SC.name = 'QNC_refractory';
        SC.param_names = {'log_QNC_Emax','log_QNC_MIC'};
        SC.prior_sd = [0.75 0.75];
        SC.prop_sd  = [0.24 0.24];

        SC.targets = struct([]);
        SC.targets(1).name         = 'QNC7D_small_n';
        SC.targets(1).regimen      = 'QNC7D';
        SC.targets(1).scenario     = 'refractory';
        SC.targets(1).x            = 3;
        SC.targets(1).n            = 3;
        SC.targets(1).response_day = 30;
        SC.targets(1).weight       = 1.0;

    case 'qnc_refractory_50_54'
        SC.name = 'qnc_refractory_50_54';
        SC.param_names = {'log_QNC_Emax','log_QNC_MIC'};
        SC.prior_sd = [0.90 0.90];
        SC.prop_sd  = [0.25 0.25];

        SC.targets = struct([]);
        SC.targets(1).name         = 'QNC_refractory_first_course';
        SC.targets(1).regimen      = 'QNC7D';
        SC.targets(1).scenario     = 'refractory';
        SC.targets(1).x            = 50;
        SC.targets(1).n            = 54;
        SC.targets(1).response_day = 30;
        SC.targets(1).weight       = 1.0;

    case 'pooled_alb_nitro_mtz5_abz5'
        SC.name = 'pooled_alb_nitro_MTZ5D_ABZ5D';
        SC.param_names = {'log_nitro_Emax','log_nitro_MIC', ...
                          'log_ABZ_Emax','log_ABZ_MIC'};
        SC.prior_sd = [0.75 0.75 0.50 0.50];
        SC.prop_sd  = [0.24 0.24 0.18 0.18];

        SC.targets = struct([]);
        SC.targets(1).name         = 'Meltzer2014_pooled_Alb_nitro';
        SC.targets(1).regimen      = 'MTZ5D_ABZ5D';
        SC.targets(1).scenario     = 'refractory';
        SC.targets(1).x            = 42;
        SC.targets(1).n            = 53;
        SC.targets(1).response_day = 30;
        SC.targets(1).weight       = 1.0;

    case 'pooled_alb_nitro_mtz7_abz7'
        SC.name = 'pooled_alb_nitro_MTZ7D_ABZ7D';
        SC.param_names = {'log_nitro_Emax','log_nitro_MIC', ...
                          'log_ABZ_Emax','log_ABZ_MIC'};
        SC.prior_sd = [0.75 0.75 0.50 0.50];
        SC.prop_sd  = [0.24 0.24 0.18 0.18];

        SC.targets = struct([]);
        SC.targets(1).name         = 'Meltzer2014_pooled_Alb_nitro';
        SC.targets(1).regimen      = 'MTZ7D_ABZ7D';
        SC.targets(1).scenario     = 'refractory';
        SC.targets(1).x            = 42;
        SC.targets(1).n            = 53;
        SC.targets(1).response_day = 30;
        SC.targets(1).weight       = 1.0;

    case 'pooled_alb_nitro_mtz14_abz7'
        SC.name = 'pooled_alb_nitro_MTZ14D_ABZ7D';
        SC.param_names = {'log_nitro_Emax','log_nitro_MIC', ...
                          'log_ABZ_Emax','log_ABZ_MIC'};
        SC.prior_sd = [0.75 0.75 0.50 0.50];
        SC.prop_sd  = [0.24 0.24 0.18 0.18];

        SC.targets = struct([]);
        SC.targets(1).name         = 'Meltzer2014_pooled_Alb_nitro';
        SC.targets(1).regimen      = 'MTZ14D_ABZ7D';
        SC.targets(1).scenario     = 'refractory';
        SC.targets(1).x            = 42;
        SC.targets(1).n            = 53;
        SC.targets(1).response_day = 30;
        SC.targets(1).weight       = 1.0;

    case 'pooled_alb_nitro_tdz2g_abz5'
        SC.name = 'pooled_alb_nitro_TDZ2G_ABZ5D';
        SC.param_names = {'log_nitro_Emax','log_nitro_MIC', ...
                          'log_ABZ_Emax','log_ABZ_MIC'};
        SC.prior_sd = [0.75 0.75 0.50 0.50];
        SC.prop_sd  = [0.24 0.24 0.18 0.18];

        SC.targets = struct([]);
        SC.targets(1).name         = 'Meltzer2014_pooled_Alb_nitro';
        SC.targets(1).regimen      = 'TDZ2G_ABZ5D';
        SC.targets(1).scenario     = 'refractory';
        SC.targets(1).x            = 42;
        SC.targets(1).n            = 53;
        SC.targets(1).response_day = 30;
        SC.targets(1).weight       = 1.0;

    case 'pooled_alb_nitro_tdz2g_abz7'
        SC.name = 'pooled_alb_nitro_TDZ2G_ABZ7D';
        SC.param_names = {'log_nitro_Emax','log_nitro_MIC', ...
                          'log_ABZ_Emax','log_ABZ_MIC'};
        SC.prior_sd = [0.75 0.75 0.50 0.50];
        SC.prop_sd  = [0.24 0.24 0.18 0.18];

        SC.targets = struct([]);
        SC.targets(1).name         = 'Meltzer2014_pooled_Alb_nitro';
        SC.targets(1).regimen      = 'TDZ2G_ABZ7D';
        SC.targets(1).scenario     = 'refractory';
        SC.targets(1).x            = 42;
        SC.targets(1).n            = 53;
        SC.targets(1).response_day = 30;
        SC.targets(1).weight       = 1.0;

    otherwise
        error('Unknown scenario_name: %s', scenario_name);
end

%% ========================================================================
%  VIRTUAL-PATIENT UNCERTAINTY SETTINGS
%  ========================================================================

% Virtual-patient uncertainty switches.
% For pure deterministic PD calibration, set these false.
% For current structured uncertainty framework, keep true.
SC.sample_parasite_block = true;
SC.sample_host_block     = true;
SC.sample_gut_block      = true;

% Standard deviations for virtual-patient uncertainty blocks.
SC.vp.parasite.rT       = 0.30;
SC.vp.parasite.K        = 0.50;
SC.vp.parasite.k_encyst = 0.30;
SC.vp.parasite.k_c      = 0.50;
SC.vp.parasite.T0       = 0.50;
SC.vp.parasite.C0       = 0.50;

SC.vp.host.k_I     = 0.50;
SC.vp.host.k_micro = 0.50;
SC.vp.host.dys     = 0.50;
SC.vp.host.rM      = 0.30;
SC.vp.host.M0_sd   = 0.20;

SC.vp.gut.k_GE      = 0.20;
SC.vp.gut.k_tr_duod = 0.20;
SC.vp.gut.V_duod_L  = 0.30;
SC.vp.gut.F_duod    = 0.35;
SC.vp.gut.alpha_sd  = 0.10;

% Success definition at response_day.
SC.success_rule = 'cleared_by_day_and_below_threshold';
% Alternative:
% SC.success_rule = 'below_threshold_only';

% Numerics / saving.
SC.dt_h = 0.25;
SC.save_outputs = true;

safe_name = regexprep(lower(SC.name), '[^a-zA-Z0-9]+', '_');
safe_name = regexprep(safe_name, '_+', '_');
safe_name = regexprep(safe_name, '^_|_$', '');

SC.output_prefix = fullfile('mcmc_outputs', ...
    sprintf('giardia_mcmc_%s_mcmc%d_vp%d', safe_name, mcmc_seed, vp_seed));

outdir = fileparts(SC.output_prefix);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% ========================================================================
%  RUN MCMC
%  ========================================================================

local_validate_scenario(SC);

p0 = giardia_params();
tvec = 0:SC.dt_h:(24*t_end_days);

% -------------------------------------------------------------------------
% Fixed virtual-patient panel.
% Holding vp_seed fixed means the same virtual population is used across
% repeated MCMC chains and across scenarios.
% -------------------------------------------------------------------------
Z = local_make_virtual_panel(n_vp, vp_seed);

% -------------------------------------------------------------------------
% MCMC proposal randomness.
% Changing mcmc_seed changes the random-walk path but not the virtual panel.
% -------------------------------------------------------------------------
rng(mcmc_seed);

if isempty(theta_init)
    theta = zeros(1, numel(SC.param_names));
else
    theta = reshape(theta_init, 1, []);
    if numel(theta) ~= numel(SC.param_names)
        error('theta_init must have length %d for this scenario.', numel(SC.param_names));
    end
end
[logpost, detail] = local_logpost(theta, p0, Z, tvec, SC);

Theta = NaN(n_iter, numel(SC.param_names));
LogPost = NaN(n_iter,1);
Accepted = false(n_iter,1);
Details = cell(n_iter,1);

for it = 1:n_iter
    theta_prop = theta + SC.prop_sd .* randn(size(theta));
    [lp_prop, detail_prop] = local_logpost(theta_prop, p0, Z, tvec, SC);

    if log(rand) < (lp_prop - logpost)
        theta = theta_prop;
        logpost = lp_prop;
        detail = detail_prop;
        Accepted(it) = true;
    end

    Theta(it,:) = theta;
    LogPost(it) = logpost;
    Details{it} = detail;

    if mod(it, max(1, floor(n_iter/10))) == 0
        fprintf('MCMC %5d/%5d | logpost %.2f | accept %.2f\n', ...
            it, n_iter, logpost, mean(Accepted(1:it)));
    end
end

%% ========================================================================
%  PACKAGE OUTPUT
%  ========================================================================

chain = struct();
chain.theta = Theta;
chain.logpost = LogPost;
chain.accepted = Accepted;
chain.accept_rate = mean(Accepted);
chain.param_names = SC.param_names;
chain.prior_sd = SC.prior_sd;
chain.prop_sd = SC.prop_sd;
chain.targets = SC.targets;
chain.scenario = SC;
chain.details = Details;
chain.n_iter = n_iter;
chain.n_vp = n_vp;
chain.t_end_days = t_end_days;
chain.mcmc_seed = mcmc_seed;
chain.vp_seed = vp_seed;
chain.when = datestr(now);

% Convenience posterior summaries after 50% burn-in.
burn = floor(n_iter/2) + 1;
S = Theta(burn:end,:);

chain.posterior_summary = table(SC.param_names(:), mean(S,1)', median(S,1)', ...
    prctile(S,2.5,1)', prctile(S,97.5,1)', exp(median(S,1))', ...
    'VariableNames', {'parameter','mean_log','median_log','q025_log','q975_log','median_multiplier'});

% Posterior-predictive target fit at posterior median.
theta_med = median(S,1);
[~, fit_detail] = local_logpost(theta_med, p0, Z, tvec, SC);
chain.posterior_median_fit = fit_detail;

disp(' ')
disp('=== Posterior summary, second half of chain ===')
disp(chain.posterior_summary)

disp(' ')
disp('=== Target fit at posterior median ===')
local_disp_target_fit(fit_detail)

if SC.save_outputs
    try
        matfile = sprintf('%s_chain.mat', SC.output_prefix);
        csvfile = sprintf('%s_posterior_summary.csv', SC.output_prefix);
        save(matfile, 'chain');
        writetable(chain.posterior_summary, csvfile);
        fprintf('\nSaved: %s\n', matfile);
        fprintf('Saved: %s\n', csvfile);
    catch ME
        warning('Could not save MCMC outputs: %s', ME.message);
    end
end

end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function local_validate_scenario(SC)
% Basic guardrails for common editing errors.

if numel(SC.param_names) ~= numel(SC.prior_sd)
    error('SC.param_names and SC.prior_sd must have the same length.');
end

if numel(SC.param_names) ~= numel(SC.prop_sd)
    error('SC.param_names and SC.prop_sd must have the same length.');
end

allowed = {'log_nitro_Emax','log_nitro_MIC', ...
           'log_ABZ_Emax','log_ABZ_MIC', ...
           'log_NTZ_Emax','log_NTZ_MIC', ...
           'log_QNC_Emax','log_QNC_MIC'};

for j = 1:numel(SC.param_names)
    if ~ismember(SC.param_names{j}, allowed)
        error('Unknown parameter name: %s', SC.param_names{j});
    end
end

required = {'name','regimen','scenario','x','n','response_day','weight'};

for k = 1:numel(SC.targets)
    for j = 1:numel(required)
        if ~isfield(SC.targets(k), required{j})
            error('SC.targets(%d) is missing field: %s', k, required{j});
        end
    end
end

end

% -------------------------------------------------------------------------
function [lp, detail] = local_logpost(theta, p0, Z, tvec, SC)
% Log posterior = log prior + weighted binomial log likelihood.

lp_prior = -0.5 * sum((theta ./ SC.prior_sd).^2) ...
           - sum(log(SC.prior_sd)) ...
           - 0.5*numel(theta)*log(2*pi);

lp_like = 0;
rows = struct([]);
row = 0;

for k = 1:numel(SC.targets)
    tar = SC.targets(k);

    if tar.weight == 0
        continue
    end

    q = local_model_response_probability(theta, p0, Z, tvec, tar, SC);
    q = min(max(q, 1e-6), 1 - 1e-6);

    % Binomial log likelihood; choose(n,x) constant omitted.
    ll = tar.x * log(q) + (tar.n - tar.x) * log(1 - q);
    lp_like = lp_like + tar.weight * ll;

    row = row + 1;
    rows(row).name = tar.name;
    rows(row).n = tar.n;
    rows(row).x = tar.x;
    rows(row).observed = tar.x / tar.n;
    rows(row).predicted = q;
    rows(row).loglike = ll;
    rows(row).weight = tar.weight;
end

lp = lp_prior + lp_like;
detail = struct('log_prior', lp_prior, ...
                'log_likelihood', lp_like, ...
                'target_rows', rows);

end

% -------------------------------------------------------------------------
function q = local_model_response_probability(theta, p0, Z, tvec, tar, SC)
% Estimate response probability from fixed virtual-patient panel.

n_vp = numel(Z.rT);
success = false(n_vp,1);

for i = 1:n_vp
    p = local_make_patient(p0, Z, i, SC);
    p = local_apply_pd_multipliers(p, theta, SC.param_names);

    opts = struct('useMicrobiome', true, 'resistant', false);

    if strcmpi(tar.scenario, 'refractory')
        opts.resistant = true;
    elseif strcmpi(tar.scenario, 'baseline')
        opts.resistant = false;
    else
        error('Unknown target scenario: %s', tar.scenario);
    end

    reg = giardia_regimens(tar.regimen, 0);
    out = giardia_simulate(tvec, reg, p, opts);

    success(i) = local_success_at_day(out, p, tar.response_day, SC.success_rule);
end

q = mean(success);

end

% -------------------------------------------------------------------------
function tf = local_success_at_day(out, p, response_day, success_rule)
% Operational clinical success proxy at a specified observation day.

Tday = interp1(out.tvec(:)/24, out.T(:), response_day, 'linear', 'extrap');

switch lower(success_rule)

    case 'cleared_by_day_and_below_threshold'
        cleared_by_day = ~isnan(out.clearance_time_h) && ...
                         out.clearance_time_h <= 24*response_day;
        tf = cleared_by_day && (Tday < p.clearance_threshold);

    case 'below_threshold_only'
        tf = Tday < p.clearance_threshold;

    otherwise
        error('Unknown success_rule: %s', success_rule);
end

end

% -------------------------------------------------------------------------
function p = local_apply_pd_multipliers(p, theta, param_names)
% Apply named log-multiplier PD parameters.

mult = struct();

for j = 1:numel(param_names)
    mult.(param_names{j}) = exp(theta(j));
end

% Nitroimidazole class: MTZ and TDZ move together.
if isfield(mult,'log_nitro_Emax')
    p.MTZ.Emax = p.MTZ.Emax * mult.log_nitro_Emax;
    p.TDZ.Emax = p.TDZ.Emax * mult.log_nitro_Emax;
end

if isfield(mult,'log_nitro_MIC')
    p.MTZ.MIC = p.MTZ.MIC * mult.log_nitro_MIC;
    p.TDZ.MIC = p.TDZ.MIC * mult.log_nitro_MIC;
end

% Albendazole block.
if isfield(mult,'log_ABZ_Emax')
    p.ABZ.Emax = p.ABZ.Emax * mult.log_ABZ_Emax;
end

if isfield(mult,'log_ABZ_MIC')
    p.ABZ.MIC = p.ABZ.MIC * mult.log_ABZ_MIC;
end

% Nitazoxanide block.
if isfield(mult,'log_NTZ_Emax')
    p.NTZ.Emax = p.NTZ.Emax * mult.log_NTZ_Emax;
end

if isfield(mult,'log_NTZ_MIC')
    p.NTZ.MIC = p.NTZ.MIC * mult.log_NTZ_MIC;
end

% Quinacrine block.
if isfield(mult,'log_QNC_Emax')
    p.QNC.Emax = p.QNC.Emax * mult.log_QNC_Emax;
end

if isfield(mult,'log_QNC_MIC')
    p.QNC.MIC = p.QNC.MIC * mult.log_QNC_MIC;
end

end

% -------------------------------------------------------------------------
function Z = local_make_virtual_panel(n_vp, seed)
% Common random numbers for virtual-patient uncertainty.

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

% Shared gut/local exposure uncertainty block.
Z.gut_GE   = randn(n_vp,1);
Z.gut_tr   = randn(n_vp,1);
Z.gut_V    = randn(n_vp,1);
Z.gut_F    = randn(n_vp,1);
Z.alpha    = randn(n_vp,1);

end

% -------------------------------------------------------------------------
function p = local_make_patient(p0, Z, i, SC)
% Build one virtual patient from fixed random numbers.

p = p0;

% A. Parasite natural history.
if SC.sample_parasite_block
    p.rT        = max(1e-6, p.rT        * exp(SC.vp.parasite.rT       * Z.rT(i)));
    p.K         = max(1e-6, p.K         * exp(SC.vp.parasite.K        * Z.K(i)));
    p.k_encyst  = max(0,    p.k_encyst  * exp(SC.vp.parasite.k_encyst * Z.k_encyst(i)));
    p.k_c_decay = max(0,    p.k_c_decay * exp(SC.vp.parasite.k_c      * Z.k_c(i)));
    p.T0        = max(1e-9, p.T0        * exp(SC.vp.parasite.T0       * Z.T0(i)));
    p.C0        = max(0,    p.C0        * exp(SC.vp.parasite.C0       * Z.C0(i)));
end

% B. Host / microbiome.
if SC.sample_host_block
    p.k_I     = max(0, p.k_I     * exp(SC.vp.host.k_I     * Z.k_I(i)));
    p.k_micro = max(0, p.k_micro * exp(SC.vp.host.k_micro * Z.k_micro(i)));
    p.dys     = max(0, p.dys     * exp(SC.vp.host.dys     * Z.dys(i)));
    p.rM      = max(0, p.rM      * exp(SC.vp.host.rM      * Z.rM(i)));
    p.M0      = min(1, max(0, p.M0 + SC.vp.host.M0_sd * Z.M0(i)));
end

% C. Shared GI/local exposure multipliers.
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
        p.exposure.alpha_plasma = min(0.8, ...
            max(0, p.exposure.alpha_plasma + SC.vp.gut.alpha_sd * Z.alpha(i)));
    end
end

end

% -------------------------------------------------------------------------
function local_disp_target_fit(detail)
% Print target rows nicely.

if ~isfield(detail, 'target_rows') || isempty(detail.target_rows)
    fprintf('No active targets.\n');
    return
end

rows = detail.target_rows;

for k = 1:numel(rows)
    fprintf('%-30s observed %.3f (%d/%d), predicted %.3f, loglike %.2f\n', ...
        rows(k).name, ...
        rows(k).observed, ...
        rows(k).x, ...
        rows(k).n, ...
        rows(k).predicted, ...
        rows(k).loglike);
end

end