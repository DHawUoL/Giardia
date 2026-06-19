function C = giardia_pk(tvec, regimen, p)
% GIARDIA_PK  Minimal oral PK with stomach, duodenum, and plasma.
%
% Output fields for each drug:
%   stom_amt   : amount in stomach (mg)
%   duod_amt   : amount in duodenum (mg)
%   plasma_amt : amount in plasma compartment (mg)
%   duod_conc  : duodenal concentration (ug/L)
%   plasma_conc: plasma concentration (ug/L)
%   eff        : currently set equal to duodenal concentration
%
% Notes:
% - Doses are added during the time-stepping loop, so repeated dosing works.
% - This is still a minimal model: no explicit dissolution, precipitation, or fed/fasted effects.

C = struct();

if ~isfield(regimen,'doses') || isempty(fieldnames(regimen.doses))
    return
end

dt = tvec(2) - tvec(1);
nT = numel(tvec);
drugs = fieldnames(regimen.doses);

for d = 1:numel(drugs)
    name = drugs{d};
    doses = regimen.doses.(name);

    A_stom   = zeros(nT,1);   % mg
    A_duod   = zeros(nT,1);   % mg
    A_plasma = zeros(nT,1);   % mg

    kGE   = p.(name).k_GE;         % stomach -> duodenum
    kAbs  = p.(name).k_abs_duod;   % duodenum -> plasma
    kTr   = p.(name).k_tr_duod;    % duodenum onward transit
    kel   = p.(name).kel;          % plasma elimination
    Vduod = p.(name).V_duod_L;     % L
    Vplas = p.(name).Vd_L;         % L
    Fduo  = p.(name).F_duod;       % fraction entering modelled oral pathway

    dose_times = [doses.time_h];
    dose_amts  = [doses.amount_mg] * Fduo;

    for i = 1:nT-1
        % Add any doses scheduled exactly at this time step
        due = abs(dose_times - tvec(i)) < dt/2;
        if any(due)
            A_stom(i) = A_stom(i) + sum(dose_amts(due));
        end

        dA_stom   = -kGE * A_stom(i);
        dA_duod   =  kGE * A_stom(i) - (kAbs + kTr) * A_duod(i);
        dA_plasma =  kAbs * A_duod(i) - kel * A_plasma(i);

        A_stom(i+1)   = max(0, A_stom(i)   + dt * dA_stom);
        A_duod(i+1)   = max(0, A_duod(i)   + dt * dA_duod);
        A_plasma(i+1) = max(0, A_plasma(i) + dt * dA_plasma);
    end

    C.(name).stom_amt    = A_stom;
    C.(name).duod_amt    = A_duod;
    C.(name).plasma_amt  = A_plasma;
    C.(name).duod_conc   = (A_duod   ./ Vduod) * 1000;  % mg/L -> ug/L
    C.(name).plasma_conc = (A_plasma ./ Vplas) * 1000;  % mg/L -> ug/L

    % Baseline effective concentration for parasite kill:
    % for now, keep it duodenum-led
    C.(name).eff = C.(name).duod_conc + p.exposure.alpha_plasma * C.(name).plasma_conc;
end

end