function plot_multichain_pairplot(chains, chain_labels, burn_frac, use_multiplier_scale)
% PLOT_MULTICHAIN_PAIRPLOT
% Pairwise posterior scatter/density diagnostic across multiple MCMC chains.
%
% Usage:
%   chains = {c1,c2,c3,c4};
%   labels = {'init 0','init +','init -','mixed init'};
%   plot_multichain_pairplot(chains, labels, 0.5, true)
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

% Check compatible chains.
for c = 2:n_chains
    if numel(chains{c}.param_names) ~= n_params || ...
       any(~strcmp(chains{c}.param_names, param_names))
        error('All chains must have the same parameter names in the same order.');
    end
end

% Pool post-burn samples.
X = [];
chain_id = [];

for c = 1:n_chains
    theta = chains{c}.theta;
    burn = floor(burn_frac * size(theta,1)) + 1;
    vals = theta(burn:end,:);

    if use_multiplier_scale
        vals = exp(vals);
    end

    X = [X; vals];
    chain_id = [chain_id; c * ones(size(vals,1),1)];
end

pretty_names = cellfun(@local_pretty_param_name, param_names, 'UniformOutput', false);

figure('Color','w','Position',[100 100 1050 950]);

% Downsample for legibility if needed.
max_points_per_chain = 1200;
keep = false(size(chain_id));

for c = 1:n_chains
    idx = find(chain_id == c);
    if numel(idx) > max_points_per_chain
        idx = idx(randperm(numel(idx), max_points_per_chain));
    end
    keep(idx) = true;
end

Xplot = X(keep,:);
chain_plot = chain_id(keep);

for i = 1:n_params
    for j = 1:n_params
        subplot(n_params, n_params, (i-1)*n_params + j)
        hold on

        if i == j
            % Diagonal: marginal densities.
            for c = 1:n_chains
                vals = X(chain_id == c, j);
                [f, xi] = ksdensity(vals);
                plot(xi, f, 'LineWidth', 1.2);
            end

            if use_multiplier_scale
                xline(1, ':k', 'HandleVisibility','off');
            else
                xline(0, ':k', 'HandleVisibility','off');
            end

            ylabel('Density')

        else
            % Off-diagonal: pairwise scatter.
            for c = 1:n_chains
                idx = chain_plot == c;
                scatter(Xplot(idx,j), Xplot(idx,i), ...
                    8, 'filled', ...
                    'MarkerFaceAlpha', 0.18, ...
                    'MarkerEdgeAlpha', 0.18);
            end

            if use_multiplier_scale
                xline(1, ':k', 'HandleVisibility','off');
                yline(1, ':k', 'HandleVisibility','off');
            else
                xline(0, ':k', 'HandleVisibility','off');
                yline(0, ':k', 'HandleVisibility','off');
            end
        end

        box on
        grid on

        if i == n_params
            xlabel(pretty_names{j}, 'Interpreter','none')
        else
            set(gca, 'XTickLabel', [])
        end

        if j == 1
            ylabel(pretty_names{i}, 'Interpreter','none')
        else
            set(gca, 'YTickLabel', [])
        end

        if i == 1 && j == n_params
            legend(chain_labels, 'Location','best')
        end
    end
end

if use_multiplier_scale
    sgtitle('Pairwise posterior structure across chains: multipliers', 'FontWeight','bold');
else
    sgtitle('Pairwise posterior structure across chains: log-multipliers', 'FontWeight','bold');
end

end

function out = local_pretty_param_name(name)

switch name
    case 'log_nitro_Emax'
        out = 'Nitro E_{max}';
    case 'log_nitro_MIC'
        out = 'Nitro conc. scale';
    case 'log_ABZ_Emax'
        out = 'ABZ E_{max}';
    case 'log_ABZ_MIC'
        out = 'ABZ conc. scale';
    case 'log_NTZ_Emax'
        out = 'NTZ E_{max}';
    case 'log_NTZ_MIC'
        out = 'NTZ conc. scale';
    case 'log_QNC_Emax'
        out = 'QNC E_{max}';
    case 'log_QNC_MIC'
        out = 'QNC conc. scale';
    otherwise
        out = strrep(name,'_',' ');
end

end

%{
chains = {c1,c2,c3,c4};
labels = {'init 0','init +','init -','mixed init'};

plot_multichain_pairplot(chains, labels, 0.5, true)
%}