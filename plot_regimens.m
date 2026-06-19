function plot_regimens()
% PLOT_REGIMENS  Quick demo: concentrations and trajectories for canonical regimens.

p = giardia_params();
tvec = 0:0.1:(24*21); % 3 weeks

names = {'MTZ5D','ABZ5D','MTZ14D_ABZ7D','TDZ2G_ONCE','QNC7D'};

figure; tiledlayout(2,1,'TileSpacing','compact');

% --- Concentration plot ---
nexttile;
hold on;
h1 = gobjects(numel(names),1);
for k = 1:numel(names)
    reg = giardia_regimens(names{k},0);
    C = giardia_pk(tvec, reg, p);

    ceff = zeros(numel(tvec),1);
    f = fieldnames(C);
    for j = 1:numel(f)
        ceff = ceff + C.(f{j}).eff;
    end

    h1(k) = plot(tvec(:)/24, ceff(:), 'DisplayName', names{k});
end
xlabel('Time (days)');
ylabel('Eff. concentration (ug/L)');
title('Effective gut-site concentration');
legend(h1, names, 'Location','northeastoutside');

% --- Trophozoite plot ---
nexttile;
hold on;
h2 = gobjects(numel(names),1);
for k = 1:numel(names)
    reg = giardia_regimens(names{k},0);
    out = giardia_simulate(tvec, reg, p, struct('useMicrobiome',true,'resistant',false));
    h2(k) = plot(tvec(:)/24, out.T(:), 'DisplayName', names{k});
end
yline(p.clearance_threshold,'k--','Clearance thr');
xlabel('Time (days)');
ylabel('Trophozoites (billions)');
title('Within-host dynamics');
legend(h2, names, 'Location','northeastoutside');

end