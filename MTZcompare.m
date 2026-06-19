function f=MTZcompare()

p = giardia_params;
p.MTZ.f_lumen = 0.01;
p.MTZ.Emax = 0.08;
p.MTZ.MIC = 20;

tvec = 0:0.1:(24*30);

reg0.name = 'NONE';
reg0.doses = struct();

reg1 = giardia_regimens('MTZ5D',0);
reg2 = giardia_regimens('MTZ14D_ABZ7D',0);

% make a long MTZ-only regimen
reg3.name = 'MTZ14D';
tM = (0:8:24*14-1);
tM = tM(tM < 24*14);
reg3.doses.MTZ = arrayfun(@(ti) struct('time_h',ti,'amount_mg',400), tM);

out0 = giardia_simulate(tvec, reg0, p, struct('useMicrobiome',true,'resistant',false));
out1 = giardia_simulate(tvec, reg1, p, struct('useMicrobiome',true,'resistant',false));
out2 = giardia_simulate(tvec, reg2, p, struct('useMicrobiome',true,'resistant',false));
out3 = giardia_simulate(tvec, reg3, p, struct('useMicrobiome',true,'resistant',false));

figure
plot(tvec/24, out0.T, 'k-', 'LineWidth', 2); hold on
plot(tvec/24, out1.T, 'b-', 'LineWidth', 2);
plot(tvec/24, out3.T, 'm-', 'LineWidth', 2);
plot(tvec/24, out2.T, 'r-', 'LineWidth', 2);
xlabel('Time (days)')
ylabel('Trophozoites (billions)')
legend('Untreated','MTZ5D','MTZ14D','MTZ14D+ABZ7D')
title('Duration vs combination')