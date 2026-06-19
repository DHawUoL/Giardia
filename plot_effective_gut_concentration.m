function plot_effective_gut_concentration(t_end_days)
% PLOT_EFFECTIVE_GUT_CONCENTRATION_STACKED
% Small-multiple plot of modelled effective gut-site concentration.
%
% One panel per regimen. Different colour per line.
% Regimen titles are placed inside the axes to save vertical space.
% One global y-axis label is used across all panels.

if nargin < 1 || isempty(t_end_days)
    t_end_days = 15;
end

p = giardia_params();
tvec = 0:0.05:(24*t_end_days);

regimen_keys = {
    'MTZ5D'
    'TDZ2G_ONCE'
    'MTZ14D_ABZ7D'
    'NTZ3D'
    'QNC7D'
};

regimen_labels = {
    'MTZ 5 days'
    'TDZ 2 g once'
    'MTZ 14 d + ABZ 7 d'
    'NTZ 3 days'
    'QNC 7 days'
};

% Distinct but still MATLAB-default-ish colours
cols = lines(numel(regimen_keys));

zero_tol = 1e-1;   % ug/L, plotting threshold only

fig = figure;
set(fig, 'Color', 'w');

tl = tiledlayout(numel(regimen_keys), 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = gobjects(numel(regimen_keys), 1);

for r = 1:numel(regimen_keys)

    reg = giardia_regimens(regimen_keys{r}, 0);
    Cpk = giardia_pk(tvec, reg, p);

    Ceff_total = zeros(numel(tvec), 1);
    drug_names = fieldnames(Cpk);

    for d = 1:numel(drug_names)
        drug = drug_names{d};
        Ceff_total = Ceff_total + Cpk.(drug).eff(:);
    end

    y = Ceff_total;
    y(y < zero_tol) = NaN;
    y = log10(y);

    ax(r) = nexttile;
    plot(tvec(:)/24, y, ...
        'LineWidth', 2, ...
        'Color', cols(r,:));

    xlim([0 t_end_days]);
    ylim([-1 5]);

    yticks(-1:2:5);
    yticklabels({'10^{-1}', '10^{1}', '10^{3}', '10^{5}'});

    grid on;
    grid minor;
    box on;

    % Regimen title inside axes, top-left
    %text(0.015, 0.82, regimen_labels{r}, ...
    text(0.015, 0.95, regimen_labels{r}, ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'FontSize', 11, ...
    'Color', 'k', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...% cols(r,:), ...
    'Margin', 3);

    if r < numel(regimen_keys)
        xticklabels([]);
    else
        xlabel('Time since treatment start (days)');
    end

    set(gca, 'FontSize', 10);
end

% Link x axes
linkaxes(ax, 'x');

% Main title
title(tl, 'Modelled effective gut-site drug exposure', ...
    'FontWeight', 'bold', ...
    'FontSize', 14);

% Shared y-label for tiledlayout if your MATLAB version supports it
try
    ylabel(tl, 'Effective gut-site concentration, C_{eff} (\mug/L)');
catch
    % Fallback for older MATLAB versions:
    han = axes(fig, 'Visible', 'off');
    han.YLabel.Visible = 'on';
    ylabel(han, 'Effective gut-site concentration, C_{eff} (\mug/L)', ...
        'FontSize', 12);
end

end