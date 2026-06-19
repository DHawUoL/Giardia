% GIARDIA_MODEL_MAIN  Starter script for Giardia within-host study.
clear; close all; clc;

addpath(pwd);

p = giardia_params();

% Visualize canonical regimens
plot_regimens();

% Monte Carlo comparison of regimens
regimens = {'MTZ5D','ABZ5D','MTZ14D_ABZ7D','TDZ2G_ONCE','NTZ3D','QNC7D'};
results = giardia_montecarlo(300, 35, regimens, 42);

% Print summary
fprintf('\n=== Monte Carlo summary (N=%d) ===\n', results.N);
for k=1:numel(regimens)
    R = results.(regimens{k});
    fprintf('%-15s  p(cleared)=%.2f  median clear (days)=%.1f\n', regimens{k}, R.p_cleared, R.median_days);
end

% Simple objective ranking
scores = zeros(numel(regimens),1);
for k=1:numel(regimens)
    scores(k) = giardia_objective(regimens{k}, p, 35, 99);
end
[ss,ix] = sort(scores);
fprintf('\nObjective (lower better) - expected time-to-clear with penalty:\n');
for j=1:numel(ix)
    fprintf('%2d) %-15s  %.2f days\n', j, regimens{ix(j)}, ss(j));
end
