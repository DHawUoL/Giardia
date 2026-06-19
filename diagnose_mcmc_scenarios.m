function D = diagnose_mcmc_scenarios(result_files, scenario_labels, burn_frac)
% DIAGNOSE_MCMC_SCENARIOS
% Basic single-chain diagnostics: acceptance rate, crude ESS, lag-1 autocorr.

if nargin < 3 || isempty(burn_frac)
    burn_frac = 0.5;
end

D = table();

for s = 1:numel(result_files)
    S = load(result_files{s}, 'chain');
    chain = S.chain;

    Xfull = chain.theta;
    burn_idx = floor(burn_frac * size(Xfull,1)) + 1;
    X = Xfull(burn_idx:end,:);

    for j = 1:numel(chain.param_names)
        x = X(:,j);

        r = table();
        r.scenario = string(scenario_labels{s});
        r.parameter = string(chain.param_names{j});
        r.accept_rate = chain.accept_rate;
        r.lag1_autocorr = local_autocorr_lag1(x);
        r.ess_approx = local_ess_initial_positive(x);
        r.n_post = numel(x);

        D = [D; r]; %#ok<AGROW>
    end
end

disp(D)
writetable(D, 'mcmc_diagnostics_summary.csv');

% Trace plots
n_scen = numel(result_files);

for s = 1:n_scen
    S = load(result_files{s}, 'chain');
    chain = S.chain;
    X = chain.theta;

    figure('Color','w','Position',[100 100 1100 700]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    for j = 1:numel(chain.param_names)
        nexttile;
        plot(X(:,j), 'LineWidth', 0.7);
        yline(0, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility','off');
        title(strrep(chain.param_names{j}, '_', '\_'));
        xlabel('Iteration');
        ylabel('Log-multiplier');
        grid on
        box on
    end

    sgtitle(sprintf('Trace plots: %s', scenario_labels{s}), ...
        'FontWeight','bold');
end

end

function r = local_autocorr_lag1(x)
x = x(:);
x = x - mean(x);
r = sum(x(1:end-1).*x(2:end)) / sum(x.^2);
end

function ess = local_ess_initial_positive(x)
% Rough ESS using initial positive autocorrelation sequence.
x = x(:);
x = x - mean(x);
n = numel(x);

maxLag = min(1000, floor(n/2));
acf = zeros(maxLag,1);

den = sum(x.^2);
for k = 1:maxLag
    acf(k) = sum(x(1:end-k).*x(1+k:end)) / den;
end

% Stop when autocorrelation becomes negative
pos = find(acf <= 0, 1, 'first');
if isempty(pos)
    K = maxLag;
else
    K = pos - 1;
end

tau = 1 + 2*sum(acf(1:K));
ess = n / max(tau, 1);
end