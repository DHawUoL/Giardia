function plot_multichain_posteriors(chains, chain_labels, burn_frac, use_multiplier_scale)
% PLOT_MULTICHAIN_POSTERIORS
% Overlay posterior density estimates from multiple MCMC chains.
%
% Usage:
%   chains = {c1,c2,c3,c4};
%   labels = {'chain 1','chain 2','chain 3','chain 4'};
%   plot_multichain_posteriors(chains, labels, 0.5, true)
%
% Inputs:
%   chains               cell array of chain structs
%   chain_labels          cell array of labels
%   burn_frac             fraction discarded as burn-in, e.g. 0.5
%   use_multiplier_scale  true: plot exp(theta); false: plot theta

if nargin < 2 || isempty(chain_labels)
    chain_labels = arrayfun(@(i) sprintf('chain %d', i), 1:numel(chains), 'UniformOutput', false);
end

if nargin < 3 || isempty(burn_frac)
    burn_frac = 0.5;
end

if nargin < 4 || isempty(use_multiplier_scale)
    use_multiplier_scale = true;
end

n_chains = numel(chains);

if n_chains == 0
    error('No chains supplied.');
end

param_names = chains{1}.param_names;
n_params = numel(param_names);

% Check chains are compatible.
for c = 2:n_chains
    if numel(chains{c}.param_names) ~= n_params || ...
       any(~strcmp(chains{c}.param_names, param_names))
        error('All chains must have the same parameter names in the same order.');
    end
end

figure('Color','w','Position',[100 100 1200 800]);

nrow = ceil(n_params/2);
ncol = min(2,n_params);

for j = 1:n_params
    subplot(nrow,ncol,j)
    hold on

    all_vals = [];

    for c = 1:n_chains
        theta = chains{c}.theta;
        burn = floor(burn_frac * size(theta,1)) + 1;
        vals = theta(burn:end,j);

        if use_multiplier_scale
            vals = exp(vals);
        end

        all_vals = [all_vals; vals(:)];

        [f, x] = ksdensity(vals);
        plot(x, f, 'LineWidth', 1.8);
    end

    if use_multiplier_scale
        xline(1, ':k', 'Baseline', ...
            'LabelVerticalAlignment','bottom', ...
            'HandleVisibility','off');
        xlabel('Multiplier')
    else
        xline(0, ':k', 'Baseline', ...
            'LabelVerticalAlignment','bottom', ...
            'HandleVisibility','off');
        xlabel('Log multiplier')
    end

    ylabel('Posterior density')
    title(local_pretty_param_name(param_names{j}), 'Interpreter','none')
    box on
    grid on

    if j == 1
        legend(chain_labels, 'Location','best')
    end

    % Sensible x limits.
    if use_multiplier_scale
        q = quantile(all_vals, [0.005 0.995]);
        if all(isfinite(q)) && q(1) < q(2)
            xlim([max(0,q(1)) q(2)])
        end
    else
        q = quantile(all_vals, [0.005 0.995]);
        if all(isfinite(q)) && q(1) < q(2)
            xlim(q)
        end
    end
end

sgtitle(sprintf('Posterior estimates across %d MCMC chains', n_chains), ...
    'FontWeight','bold');

end

function out = local_pretty_param_name(name)

switch name
    case 'log_nitro_Emax'
        out = 'Nitroimidazole E_{max}';
    case 'log_nitro_MIC'
        out = 'Nitroimidazole concentration scale';
    case 'log_ABZ_Emax'
        out = 'Albendazole E_{max}';
    case 'log_ABZ_MIC'
        out = 'Albendazole concentration scale';
    case 'log_NTZ_Emax'
        out = 'Nitazoxanide E_{max}';
    case 'log_NTZ_MIC'
        out = 'Nitazoxanide concentration scale';
    case 'log_QNC_Emax'
        out = 'Quinacrine E_{max}';
    case 'log_QNC_MIC'
        out = 'Quinacrine concentration scale';
    otherwise
        out = strrep(name,'_',' ');
end

end
%{
chains = {c1,c2,c3,c4};
labels = {'init 0','init +','init -','mixed init'};

plot_multichain_posteriors(chains, labels, 0.5, true)
%}