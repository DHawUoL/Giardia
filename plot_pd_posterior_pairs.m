function plot_pd_posterior_pairs(chain, burn_frac)
% PLOT_PD_POSTERIOR_PAIRS
% Pairwise posterior plots for PD log-multipliers / multipliers.
%
% Example:
%   plot_pd_posterior_pairs(chain, 0.5)

if nargin < 2 || isempty(burn_frac)
    burn_frac = 0.5;
end

% Extract samples
if isfield(chain, 'theta')
    samples = chain.theta;
elseif isfield(chain, 'samples')
    samples = chain.samples;
else
    error('Could not find chain.theta or chain.samples.');
end

if isfield(chain, 'param_names')
    names = chain.param_names;
else
    names = {'log_nitro_Emax','log_nitro_MIC','log_ABZ_Emax','log_ABZ_MIC'};
end

n_iter = size(samples,1);
burn = floor(burn_frac*n_iter) + 1;
Slog = samples(burn:end,:);
S = exp(Slog);  % multiplier scale

% Find columns
iNE = find(strcmp(names, 'log_nitro_Emax'));
iNM = find(strcmp(names, 'log_nitro_MIC'));
iAE = find(strcmp(names, 'log_ABZ_Emax'));
iAM = find(strcmp(names, 'log_ABZ_MIC'));

% Helpful pair list
pairs = [
    iNE iAE
    iNE iNM
    iAE iAM
    iNE iAM
];

pair_labels = {
    'Nitro E_{max} vs ABZ E_{max}'
    'Nitro E_{max} vs Nitro MIC'
    'ABZ E_{max} vs ABZ MIC'
    'Nitro E_{max} vs ABZ MIC'
};

figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for k = 1:size(pairs,1)
    nexttile;
    x = S(:,pairs(k,1));
    y = S(:,pairs(k,2));

    scatter(x, y, 12, 'filled', ...
        'MarkerFaceAlpha', 0.18, ...
        'MarkerEdgeAlpha', 0.18);
    hold on;

    xline(1, ':k', 'Baseline');
    yline(1, ':k', 'Baseline');

    r = corr(log(x), log(y), 'Rows', 'complete');

    xlabel(pretty_name(names{pairs(k,1)}));
    ylabel(pretty_name(names{pairs(k,2)}));
    title(sprintf('%s, corr(log)=%.2f', pair_labels{k}, r), ...
        'Interpreter','tex');

    box on;
    grid on;
end

sgtitle('Pairwise posterior structure: PD multipliers');

end

% -------------------------------------------------------------------------
function out = pretty_name(name)

switch name
    case 'log_nitro_Emax'
        out = 'Nitro E_{max} multiplier';
    case 'log_nitro_MIC'
        out = 'Nitro MIC multiplier';
    case 'log_ABZ_Emax'
        out = 'ABZ E_{max} multiplier';
    case 'log_ABZ_MIC'
        out = 'ABZ MIC multiplier';
    otherwise
        out = strrep(name, '_', '\_');
end

end