function T = summarise_pd_sensitivity(result_files, scenario_labels, burn_frac)
% PLOT_PD_MULTIPLIER_BOXPLOTS
% Posterior multiplier distributions by assumed Alb + nitroimidazole regimen.
%
% Shows posterior samples directly, rather than intervals around medians.

if nargin < 3 || isempty(burn_frac)
    burn_frac = 0.5;
end

Tall = table();

for s = 1:numel(result_files)
    S = load(result_files{s}, 'chain');
    chain = S.chain;

    X = chain.theta;
    burn_idx = floor(burn_frac * size(X,1)) + 1;
    X = X(burn_idx:end,:);

    for j = 1:numel(chain.param_names)
        m = exp(X(:,j));

        tmp = table();
        tmp.scenario = repmat(string(scenario_labels{s}), numel(m), 1);
        tmp.parameter = repmat(string(local_pretty_name(chain.param_names{j})), numel(m), 1);
        tmp.multiplier = m;

        Tall = [Tall; tmp]; %#ok<AGROW>
    end
end

T = Tall;
writetable(T, 'pd_multiplier_posterior_samples_long.csv');

params = unique(T.parameter, 'stable');

figure('Color','w','Position',[100 100 1300 820]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for p = 1:numel(params)
    nexttile;
    idx = T.parameter == params(p);
    Tp = T(idx,:);

    boxchart(categorical(Tp.scenario), Tp.multiplier);
    yline(1, ':', 'Color', [0.35 0.35 0.35], ...
        'LineWidth', 1.3, 'HandleVisibility','off');

    title(params(p), 'Interpreter','tex');
    ylabel('Posterior multiplier');
    xlabel('');
    grid on
    box on

    xtickangle(25)
end

sgtitle('Posterior multiplier distributions by assumed Alb + nitroimidazole regimen', ...
    'FontWeight','bold');

end

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