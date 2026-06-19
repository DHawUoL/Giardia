function plot_pd_schematic(export_name)
% PLOT_PD_SCHEMATIC
% Illustrative schematic for Emax and MIC-like concentration scale.
%
% Example:
%   plot_pd_schematic('pd_schematic.png')

if nargin < 1 || isempty(export_name)
    export_name = 'pd_schematic.png';
end

C = logspace(-2, 1, 400);
h = 1.5;

fig = figure('Color','w','Position',[100 100 600 520]);%[100 100 950 520]
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ------------------------------------------------------------
% Panel 1: Emax
% ------------------------------------------------------------
nexttile;
hold on;

MIC0 = 1.0;
Evals = [0.6 1.0 1.5];

for Emax = Evals
    y = Emax .* (C.^h) ./ (MIC0^h + C.^h);
    semilogx(C, y, 'LineWidth', 2);
end

xline(MIC0, ':k', 'MIC', ...
    'LabelVerticalAlignment','bottom', ...
    'LabelOrientation','horizontal');

xlabel('Effective concentration, C');
ylabel('Drug kill rate');
title('E_{max}: maximum achievable drug effect');
legend({'Lower E_{max}','Baseline E_{max}','Higher E_{max}'}, ...
    'Location','southeast');

text(0.03, 0.90, ...
    {'Changing E_{max}', 'raises or lowers the plateau'}, ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'BackgroundColor','w', ...
    'EdgeColor', 'k', ...
    'Margin',4);

box on; grid on;
set(gca,'FontSize',12);

% ------------------------------------------------------------
% Panel 2: MIC / concentration scale
% ------------------------------------------------------------
nexttile;
hold on;

Emax0 = 1.2;
MICvals = [0.5 1.0 2.0];

for MIC = MICvals
    y = Emax0 .* (C.^h) ./ (MIC^h + C.^h);
    semilogx(C, y, 'LineWidth', 2);
end

xline(1.0, ':k', 'Baseline MIC', ...
    'LabelVerticalAlignment','bottom', ...
    'LabelOrientation','horizontal');

xlabel('Effective concentration, C');
ylabel('Drug kill rate');
title('MIC parameter: concentration scale for effect');
legend({'Lower MIC','Baseline MIC','Higher MIC'}, ...
    'Location','southeast');

text(0.03, 0.90, ...
    {'Changing MIC shifts the curve', 'left or right along concentration'}, ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'BackgroundColor','w', ...
    'EdgeColor', 'k', ...
    'Margin',4);

box on; grid on;
set(gca,'FontSize',12);

exportgraphics(fig, export_name, 'Resolution', 300);
fprintf('Saved: %s\n', export_name);

end