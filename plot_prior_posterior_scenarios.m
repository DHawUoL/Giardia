function plot_prior_posterior_scenarios(result_files, scenario_labels)%, prior_sd, burn_frac)
% PLOT_PRIOR_POSTERIOR_SCENARIOS
% 2x2 panel: prior (black dashed) + posterior density for each scenario.
%{
% Example:
 result_files = {
   'mcmc_sensitivity_outputs/pooled_alb_nitro_mtz5_abz5_chain.mat'
   'mcmc_sensitivity_outputs/pooled_alb_nitro_mtz7_abz7_chain.mat'
   'mcmc_sensitivity_outputs/pooled_alb_nitro_mtz14_abz7_chain.mat'
   'mcmc_sensitivity_outputs/pooled_alb_nitro_tdz2g_abz5_chain.mat'
   'mcmc_sensitivity_outputs/pooled_alb_nitro_tdz2g_abz7_chain.mat'
   'qnc_refractory_50_54_chain'
   };

 scenario_labels = {
   'MTZ5D + ABZ5D'
   'MTZ7D + ABZ7D'
   'MTZ14D + ABZ7D'
   'TDZ2G + ABZ5D'
   'TDZ2G + ABZ7D'
   'QNC7D'

   };
%}
% prior_sd = [0.75 0.75 0.50 0.50];
% plot_prior_posterior_scenarios(result_files, scenario_labels, prior_sd, 0.5)

if nargin < 3 || isempty(prior_sd)
    prior_sd = [0.75 0.75 0.50 0.50];
end
if nargin < 4 || isempty(burn_frac)
    burn_frac = 0.5;
end

n_scen = numel(result_files);
if numel(scenario_labels) ~= n_scen
    error('scenario_labels must have same length as result_files.');
end

% Load first chain to get parameter names
tmp = load(result_files{1}, 'chain');
chain0 = tmp.chain;

if isfield(chain0, 'samples')
    samples0 = chain0.samples;
elseif isfield(chain0, 'theta')
    samples0 = chain0.theta;
elseif isfield(chain0, 'chain')
    samples0 = chain0.chain;
else
    error('Could not find MCMC samples in chain struct.');
end

if isfield(chain0, 'param_names')
    param_names = chain0.param_names;
else
    param_names = {'log_nitro_Emax','log_nitro_MIC', ...
                   'log_ABZ_Emax','log_ABZ_MIC'};
end

n_param = size(samples0,2);
if n_param ~= 4
    error('This function currently expects 4 calibrated parameters.');
end
if numel(prior_sd) ~= 4
    error('prior_sd must have length 4.');
end

% Load all posterior samples
post_samples = cell(n_scen,1);
all_vals = cell(n_param,1);

for s = 1:n_scen
    S = load(result_files{s}, 'chain');
    ch = S.chain;

    if isfield(ch, 'samples')
        X = ch.samples;
    elseif isfield(ch, 'theta')
        X = ch.theta;
    elseif isfield(ch, 'chain')
        X = ch.chain;
    else
        error('Could not find MCMC samples in chain struct for scenario %d.', s);
    end

    burn_idx = max(1, floor(burn_frac*size(X,1)) + 1);
    Xpost = X(burn_idx:end, :);

    post_samples{s} = Xpost;

    for j = 1:n_param
        all_vals{j} = [all_vals{j}; Xpost(:,j)];
    end
end

% Colours: MATLAB lines palette
cols = lines(n_scen);

figure('Color','w','Position',[100 100 1200 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for j = 1:n_param
    nexttile;
    hold on;

    % Plot range
    pooled = all_vals{j};
    lo = min([pooled; -3.5*prior_sd(j)]) - 0.2;
    hi = max([pooled;  3.5*prior_sd(j)]) + 0.2;
    xx = linspace(lo, hi, 500);

    % Prior density
    prior_y = normpdf(xx, 0, prior_sd(j));
    plot(xx, prior_y, 'k--', 'LineWidth', 2.2, 'DisplayName', 'Prior');

    % Posteriors
    for s = 1:n_scen
        Xpost = post_samples{s};
        [f, xgrid] = ksdensity(Xpost(:,j), xx);
        plot(xgrid, f, 'LineWidth', 2, 'Color', cols(s,:), ...
            'DisplayName', scenario_labels{s});
    end

    % Baseline reference
    xline(0, ':', 'Color', [0.4 0.4 0.4], ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

    title(local_pretty_name(param_names{j}), 'Interpreter', 'tex');
    xlabel('Log-multiplier');
    ylabel('Density');
    box on
    grid on
end

lgd = legend('Location','southoutside', 'NumColumns', min(3, n_scen+1));
lgd.Layout.Tile = 'south';

sgtitle('Prior and posterior distributions across pooled Alb + nitroimidazole scenarios', ...
    'FontWeight','bold');

end

% -------------------------------------------------------------------------
function out = local_pretty_name(name)

switch lower(name)
    case 'log_nitro_emax'
        out = 'Nitroimidazole E_{max}';
    case 'log_nitro_mic'
        out = 'Nitroimidazole concentration scale';
    case 'log_abz_emax'
        out = 'Albendazole E_{max}';
    case 'log_abz_mic'
        out = 'Albendazole concentration scale';
    otherwise
        out = name;
end

end

%{
% How to call:

result_files = {
    'mcmc_sensitivity_outputs/pooled_alb_nitro_mtz5_abz5_chain.mat'
    'mcmc_sensitivity_outputs/pooled_alb_nitro_mtz7_abz7_chain.mat'
    'mcmc_sensitivity_outputs/pooled_alb_nitro_mtz14_abz7_chain.mat'
    'mcmc_sensitivity_outputs/pooled_alb_nitro_tdz2g_abz5_chain.mat'
    'mcmc_sensitivity_outputs/pooled_alb_nitro_tdz2g_abz7_chain.mat'
    };

scenario_labels = {
    'MTZ5D + ABZ5D'
    'MTZ7D + ABZ7D'
    'MTZ14D + ABZ7D'
    'TDZ2G + ABZ5D'
    'TDZ2G + ABZ7D'
    };

prior_sd = [0.75 0.75 0.50 0.50];

plot_prior_posterior_scenarios(result_files, scenario_labels, prior_sd, 0.5);
%}