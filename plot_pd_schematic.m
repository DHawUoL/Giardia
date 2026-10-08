function plot_pd_schematic(export_name)
% PLOT_PD_SCHEMATIC
% Illustrative schematic for Emax and concentration-scale parameters.
%
% 3-panel layout:
%   Top (spanning both columns): C_eff(t) snippet and resulting drug-kill
%                                curves for two PD parameter settings
%   Bottom-left:                 Emax schematic
%   Bottom-right:                concentration-scale schematic
%
% Example:
%   plot_pd_schematic('pd_schematic.png')

if nargin < 1 || isempty(export_name)
    export_name = 'pd_schematic.png';
end

%% ------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------
font_size = 14;
lw = 3;

Cmap = lines(7);
col_ceff = [0.15 0.15 0.15];
col_1    = Cmap(1,:);
col_2    = Cmap(7,:);

h = 1.5;
C = logspace(-2, 1, 400);

fig = figure('Color','w','Position',[100 100 1000 720]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

%% ------------------------------------------------------------------------
% Panel 1: dynamic example
% -------------------------------------------------------------------------
ax1 = nexttile([1 2]);
hold(ax1,'on');

t = linspace(0, 18, 500);  % hours

% Stylised exposure snippet
Ceff = 0.08 ...
     + 2.4*exp(-0.30*(t-1.5).^2) ...
     + 0.9*exp(-0.10*(t-7.5).^2);

% Two illustrative PD settings
Emax1 = 1.0;  MIC1 = 0.8;
Emax2 = 1.6;  MIC2 = 1.6;

kdrug1 = Emax1 .* (Ceff.^h) ./ (MIC1^h + Ceff.^h);
kdrug2 = Emax2 .* (Ceff.^h) ./ (MIC2^h + Ceff.^h);

yyaxis(ax1,'left')
h1 = plot(ax1, t, Ceff, '--', 'Color', col_ceff, 'LineWidth', lw);
ylabel(ax1, 'Effective concentration, C_{eff}(t)', 'FontSize', font_size);
set(ax1,'YColor',col_ceff)

yyaxis(ax1,'right')
h2 = plot(ax1, t, kdrug1, '-', 'Color', col_1, 'LineWidth', lw);
h3 = plot(ax1, t, kdrug2, '-', 'Color', col_2, 'LineWidth', lw);
ylabel(ax1, 'Drug kill rate, k_{drug}(t)', 'FontSize', font_size);
set(ax1,'YColor',[0.2 0.2 0.2])

xlabel(ax1, 'Time (hours)', 'FontSize', font_size);
title(ax1, 'Illustrative pharmacodynamic response to a common exposure profile', ...
    'FontSize', font_size, 'FontWeight', 'bold');

legend(ax1, [h1 h2 h3], { ...
    'C_{eff}(t)', ...
    sprintf('k_{drug}(t): E_{max} = %.1f, concentration scale = %.1f', Emax1, MIC1), ...
    sprintf('k_{drug}(t): E_{max} = %.1f, concentration scale = %.1f', Emax2, MIC2)}, ...
    'Location','northeast', 'FontSize', font_size-1);

box(ax1,'on'); 
grid(ax1,'on');
set(ax1,'FontSize',font_size);

%% ------------------------------------------------------------------------
% Panel 2: Emax schematic
% -------------------------------------------------------------------------
ax2 = nexttile;
hold(ax2,'on');

MIC0 = 1.0;
Evals = [0.6 1.0 1.5];

hE = gobjects(numel(Evals),1);
for i = 1:numel(Evals)
    Emax = Evals(i);
    y = Emax .* (C.^h) ./ (MIC0^h + C.^h);
    hE(i) = semilogx(ax2, C, y, 'LineWidth', lw);
end

xlabel(ax2, 'Effective concentration, C', 'FontSize', font_size);
ylabel(ax2, 'Drug kill rate', 'FontSize', font_size);
title(ax2, 'Effect of E_{max}', ...
    'FontSize', font_size, 'FontWeight', 'bold');

legend(ax2, hE, ...
    {sprintf('E_{max} = %.1f', Evals(1)), ...
     sprintf('E_{max} = %.1f', Evals(2)), ...
     sprintf('E_{max} = %.1f', Evals(3))}, ...
    'Location','southeast', 'FontSize', font_size-1);

box(ax2,'on'); 
grid(ax2,'on');
set(ax2,'FontSize',font_size);

%% ------------------------------------------------------------------------
% Panel 3: concentration-scale schematic
% -------------------------------------------------------------------------
ax3 = nexttile;
hold(ax3,'on');

Emax0 = 1.2;
MICvals = [0.5 1.0 2.0];

hM = gobjects(numel(MICvals),1);
for i = 1:numel(MICvals)
    MIC = MICvals(i);
    y = Emax0 .* (C.^h) ./ (MIC^h + C.^h);
    hM(i) = semilogx(ax3, C, y, 'LineWidth', lw);
end

xlabel(ax3, 'Effective concentration, C', 'FontSize', font_size);
ylabel(ax3, 'Drug kill rate', 'FontSize', font_size);
title(ax3, 'Effect of concentration scale', ...
    'FontSize', font_size, 'FontWeight', 'bold');

legend(ax3, hM, ...
    {sprintf('Scale = %.1f', MICvals(1)), ...
     sprintf('Scale = %.1f', MICvals(2)), ...
     sprintf('Scale = %.1f', MICvals(3))}, ...
    'Location','southeast', 'FontSize', font_size-1);

box(ax3,'on'); 
grid(ax3,'on');
set(ax3,'FontSize',font_size);

%% ------------------------------------------------------------------------
% Export
% -------------------------------------------------------------------------
exportgraphics(fig, export_name, 'Resolution', 300);
fprintf('Saved: %s\n', export_name);

end