function plot_effective_gut_concentration_main(t_end_days)
% PLOT_EFFECTIVE_GUT_CONCENTRATION_MAIN
% Plot effective gut-site concentration only for the main study regimens:
%   (1) MTZ14D + ABZ7D
%   (2) QNC7D
%
% Left panel: MTZ14D + ABZ7D, with separate lines for MTZ and ABZ
% Right panel: QNC7D
%
% Effective concentration is shown on a log10 scale for readability.

if nargin < 1 || isempty(t_end_days)
    t_end_days = 15;
end

%% ------------------------------------------------------------------------
% Settings
% -------------------------------------------------------------------------
font_size = 14;
lw = 3;

p = giardia_params();
tvec = 0:0.05:(24*t_end_days);
t_days = tvec(:) / 24;

C = lines(7);
col_mtz = C(1,:);
col_abz = C(2,:);
col_qnc = C(7,:);

zero_tol = 1e-1;   % plotting threshold only (ug/L)

%% ------------------------------------------------------------------------
% Simulate PK outputs
% -------------------------------------------------------------------------
reg_combo = giardia_regimens('MTZ14D_ABZ7D', 0);
Cpk_combo = giardia_pk(tvec, reg_combo, p);

reg_qnc = giardia_regimens('QNC7D', 0);
Cpk_qnc = giardia_pk(tvec, reg_qnc, p);

% Extract effective concentrations
Ceff_mtz = Cpk_combo.MTZ.eff(:);
Ceff_abz = Cpk_combo.ABZ.eff(:);
Ceff_qnc = Cpk_qnc.QNC.eff(:);

% Apply plotting threshold and log-transform for display
y_mtz = Ceff_mtz;
y_abz = Ceff_abz;
y_qnc = Ceff_qnc;

y_mtz(y_mtz < zero_tol) = NaN;
y_abz(y_abz < zero_tol) = NaN;
y_qnc(y_qnc < zero_tol) = NaN;

y_mtz = log10(y_mtz);
y_abz = log10(y_abz);
y_qnc = log10(y_qnc);

%% ------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
fig = figure('Color','w','Position',[100 100 1100 420]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ---------------- Left panel: MTZ + ABZ ----------------
ax1 = nexttile;
hold(ax1,'on');

h_mtz = plot(ax1, t_days, y_mtz, '-', ...
    'LineWidth', lw, ...
    'Color', col_mtz);

h_abz = plot(ax1, t_days, y_abz, '-', ...
    'LineWidth', lw, ...
    'Color', col_abz);

xlim(ax1, [0 t_end_days]);
ylim(ax1, [-1 5]);

yticks(ax1, -1:2:5);
yticklabels(ax1, {'10^{-1}','10^{1}','10^{3}','10^{5}'});

xlabel(ax1, 'Time since treatment start (days)', 'FontSize', font_size);
title(ax1, 'MTZ 14 d + ABZ 7 d', 'FontSize', font_size, 'FontWeight', 'bold');

legend(ax1, [h_mtz h_abz], {'Metronidazole','Albendazole'}, ...
    'Location','northeast', ...
    'FontSize', font_size-1);

grid(ax1,'on');
grid(ax1,'minor');
box(ax1,'on');
set(ax1,'FontSize',font_size);

% End-of-dosing markers: excluded from legend
xline(ax1, 7, ':', ...
    'Color', col_abz, ...
    'LineWidth', 1.4, ...
    'HandleVisibility','off');

xline(ax1, 14, ':', ...
    'Color', col_mtz, ...
    'LineWidth', 1.4, ...
    'HandleVisibility','off');

% Black readable labels, separate from the xline objects
yl = ylim(ax1);
y_label = yl(1) + 0.18*(yl(2)-yl(1));

text(ax1, 7, y_label, 'ABZ ends', ...
    'Color','k', ...
    'FontSize', font_size-2, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Rotation',90, ...
    'BackgroundColor','none', ...
    'Margin',2, ...
    'Clipping','on');

text(ax1, 14, y_label, 'MTZ ends', ...
    'Color','k', ...
    'FontSize', font_size-2, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Rotation',90, ...
    'BackgroundColor','none', ...
    'Margin',2, ...
    'Clipping','on');

% ---------------- Right panel: QNC ----------------
ax2 = nexttile;
hold(ax2,'on');

h_qnc = plot(ax2, t_days, y_qnc, '-', ...
    'LineWidth', lw, ...
    'Color', col_qnc);

xlim(ax2, [0 t_end_days]);
ylim(ax2, [-1 5]);

yticks(ax2, -1:2:5);
yticklabels(ax2, {'10^{-1}','10^{1}','10^{3}','10^{5}'});

xlabel(ax2, 'Time since treatment start (days)', 'FontSize', font_size);
title(ax2, 'QNC 7 d', 'FontSize', font_size, 'FontWeight', 'bold');

legend(ax2, h_qnc, {'Quinacrine'}, ...
    'Location','northeast', ...
    'FontSize', font_size-1);

grid(ax2,'on');
grid(ax2,'minor');
box(ax2,'on');
set(ax2,'FontSize',font_size);

xline(ax2, 7, ':', ...
    'Color', col_qnc, ...
    'LineWidth', 1.4, ...
    'HandleVisibility','off');

yl = ylim(ax2);
y_label = yl(1) + 0.18*(yl(2)-yl(1));

text(ax2, 7, y_label, 'QNC ends', ...
    'Color','k', ...
    'FontSize', font_size-2, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Rotation',90, ...
    'BackgroundColor','none', ...
    'Margin',2, ...
    'Clipping','on');

% Shared y-label
try
    ylabel(tl, 'Effective gut-site concentration, C_{eff} (\mug/L)', ...
        'FontSize', font_size);
catch
    han = axes(fig, 'Visible', 'off');
    han.YLabel.Visible = 'on';
    ylabel(han, 'Effective gut-site concentration, C_{eff} (\mug/L)', ...
        'FontSize', font_size);
end

end