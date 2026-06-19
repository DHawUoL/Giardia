function plot_prior_posterior(chain, prior_sd, burn_frac, use_multiplier_scale)
% PLOT_PRIOR_POSTERIOR
% Plot prior vs posterior for fitted MCMC parameters.
%
% Inputs
%   chain                 output struct from giardia_montecarlo(...)
%   prior_sd              vector of prior SDs on log-multipliers
%   burn_frac             fraction of chain to discard as burn-in (default 0.5)
%   use_multiplier_scale  if true, plot exp(theta); else plot theta
%
% Example:
%   chain = giardia_montecarlo(500,50,35,42);
%   prior_sd = [0.75 0.75 0.50 0.50];
%   plot_prior_posterior(chain, prior_sd, 0.5, true)

if nargin < 2 || isempty(prior_sd)
    prior_sd = [0.75 0.75 0.50 0.50];
end
if nargin < 3 || isempty(burn_frac)
    burn_frac = 0.5;
end
if nargin < 4 || isempty(use_multiplier_scale)
    use_multiplier_scale = false;
end

% ---- extract samples robustly ----
if isfield(chain, 'samples')
    samples = chain.samples;
elseif isfield(chain, 'theta')
    samples = chain.theta;
elseif isfield(chain, 'chain')
    samples = chain.chain;
else
    error('Could not find MCMC samples in chain struct.');
end

[n_iter, n_param] = size(samples);

if isfield(chain, 'param_names')
    param_names = chain.param_names;
else
    % fallback names
    default_names = {'log_nitro_Emax','log_nitro_MIC', ...
                     'log_ABZ_Emax','log_ABZ_MIC', ...
                     'log_NTZ_Emax','log_NTZ_MIC', ...
                     'log_QNC_Emax','log_QNC_MIC'};
    param_names = default_names(1:n_param);
end

if numel(prior_sd) ~= n_param
    error('Length of prior_sd (%d) must match number of parameters (%d).', ...
          numel(prior_sd), n_param);
end

% ---- burn-in ----
burn_idx = max(1, floor(burn_frac * n_iter) + 1);
post = samples(burn_idx:end, :);

% ---- plotting layout ----
ncol = 2;
nrow = ceil(n_param / ncol);

figure;
tiledlayout(nrow, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:n_param
    nexttile;
    hold on;

    post_j = post(:, j);

    if use_multiplier_scale
        % Posterior on multiplier scale
        x_post = exp(post_j);

        % Prior on multiplier scale is lognormal with meanlog=0, sdlog=prior_sd(j)
        lo = max(0, prctile(x_post, 1) * 0.7);
        hi = prctile(x_post, 99) * 1.3;
        if hi <= 0 || ~isfinite(hi)
            hi = 3;
        end
        if lo <= 0
            lo = 1e-3;
        end

        xx = linspace(lo, hi, 400);
        prior_y = lognpdf(xx, 0, prior_sd(j));

        % posterior density
        [f_post, x_grid] = ksdensity(x_post, xx, 'Support', 'positive');

        plot(xx, prior_y, '--', 'LineWidth', 1.8);
        plot(x_grid, f_post, '-', 'LineWidth', 2.2);

        xline(1, ':k', 'Baseline', 'LabelVerticalAlignment','bottom');

        xlabel('Multiplier');
        title(local_pretty_name(param_names{j}), 'Interpreter', 'none');

    else
        % Posterior on log scale
        lo = min(prctile(post_j,1), -3*prior_sd(j)) - 0.2;
        hi = max(prctile(post_j,99),  3*prior_sd(j)) + 0.2;
        xx = linspace(lo, hi, 400);

        prior_y = normpdf(xx, 0, prior_sd(j));
        [f_post, x_grid] = ksdensity(post_j, xx);

        plot(xx, prior_y, '--', 'LineWidth', 1.8);
        plot(x_grid, f_post, '-', 'LineWidth', 2.2);

        xline(0, ':k', 'Baseline', 'LabelVerticalAlignment','bottom');

        xlabel('log-multiplier');
        title(local_pretty_name(param_names{j}), 'Interpreter', 'none');
    end

    ylabel('Density');
    box on;
    grid on;
    legend({'Prior','Posterior'}, 'Location', 'best');
end

if use_multiplier_scale
    sgtitle('Prior vs posterior: parameter multipliers');
else
    sgtitle('Prior vs posterior: log-multipliers');
end

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
    case 'log_ntz_emax'
        out = 'Nitazoxanide E_{max}';
    case 'log_ntz_mic'
        out = 'Nitazoxanide concentration scale';
    case 'log_qnc_emax'
    out = 'Quinacrine E_{max}';
case 'log_qnc_mic'
    out = 'Quinacrine concentration scale';
    otherwise
        out = name;
end

end