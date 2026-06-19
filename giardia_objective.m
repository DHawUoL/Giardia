function score = giardia_objective(regimen_name, p, t_end_days, seed)
% GIARDIA_OBJECTIVE  Simple objective: expected time-to-clear (days), penalize non-clearance.
if nargin<3, t_end_days=28; end
if nargin<4, seed=1; end
rng(seed);
tvec = 0:0.25:(24*t_end_days);
res = zeros(200,1);
for i=1:200
    sim_p = p;
    % variability
    sim_p.rT       = p.rT * exp(0.3*randn);
    sim_p.K        = p.K  * exp(0.5*randn);
    resist = rand < p.resistance.frac;
    reg = giardia_regimens(regimen_name,0);
    out = giardia_simulate(tvec, reg, sim_p, struct('useMicrobiome',true,'resistant',resist));
    if isnan(out.clearance_time_h)
        res(i) = t_end_days + 7; % penalty
    else
        res(i) = out.clearance_time_h/24;
    end
end
score = mean(res);
end
