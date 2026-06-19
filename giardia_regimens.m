function regimen = giardia_regimens(which, start_h)
% GIARDIA_REGIMENS  Build canonical Giardia dosing schedules.
%
% Usage:
%   regimen = giardia_regimens('MTZ5D', 0);
%   regimen = giardia_regimens('MTZ7D_ABZ7D', 0);
%   regimen = giardia_regimens('TDZ2G_ABZ7D', 0);
%
% Naming convention:
%   MTZ = metronidazole
%   TDZ = tinidazole
%   ABZ = albendazole
%   NTZ = nitazoxanide
%   QNC = quinacrine
%
% Doses are deliberately rough clinical proxies for modelling.
% Exact study regimens should be checked before treating any schedule as
% a direct reproduction of a published protocol.

if nargin < 2 || isempty(start_h)
    start_h = 0;
end

% Normalise names so MTZ14D+ABZ7D, MTZ14D_ABZ7D, mtz14d-abz7d all work.
key = upper(which);
key = regexprep(key, '[^A-Z0-9]+', '_');
key = regexprep(key, '_$', '');

regimen = struct();
regimen.name = key;
regimen.doses = struct();

switch key

    %% --------------------------------------------------------------------
    %  Single-drug reference regimens
    %  --------------------------------------------------------------------

    case 'MTZ5D'
        % Metronidazole 400 mg TID for 5 days.
        regimen.doses.MTZ = local_doses(start_h, 8, 5, 400);

    case 'MTZ7D'
        % Metronidazole 400 mg TID for 7 days.
        regimen.doses.MTZ = local_doses(start_h, 8, 7, 400);

    case 'MTZ10D'
        % Metronidazole 400 mg TID for 10 days.
        regimen.doses.MTZ = local_doses(start_h, 8, 10, 400);

    case 'MTZ14D'
        % Metronidazole 400 mg TID for 14 days.
        regimen.doses.MTZ = local_doses(start_h, 8, 14, 400);

    case 'ABZ5D'
        % Albendazole 400 mg once daily for 5 days.
        regimen.doses.ABZ = local_doses(start_h, 24, 5, 400);

    case 'ABZ7D'
        % Albendazole 400 mg once daily for 7 days.
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);

    case 'ABZ10D'
        % Albendazole 400 mg once daily for 10 days.
        regimen.doses.ABZ = local_doses(start_h, 24, 10, 400);

    case 'ABZ14D'
        % Albendazole 400 mg once daily for 14 days.
        regimen.doses.ABZ = local_doses(start_h, 24, 14, 400);

    case {'TDZ2G', 'TDZ2G_ONCE', 'TDZ2G_SINGLE'}
        % Tinidazole single 2 g dose.
        regimen.doses.TDZ = local_doses(start_h, 24, 1, 2000);
        regimen.doses.TDZ = regimen.doses.TDZ(1);

    case 'NTZ3D'
        % Nitazoxanide 500 mg BID for 3 days.
        regimen.doses.NTZ = local_doses(start_h, 12, 3, 500);

    case 'NTZ7D'
        % Nitazoxanide 500 mg BID for 7 days.
        regimen.doses.NTZ = local_doses(start_h, 12, 7, 500);

    case 'QNC5D'
        % Quinacrine 100 mg TID for 5 days.
        regimen.doses.QNC = local_doses(start_h, 8, 5, 100);

    case 'QNC7D'
        % Quinacrine 100 mg TID for 7 days.
        regimen.doses.QNC = local_doses(start_h, 8, 7, 100);

    case 'QNC10D'
        % Quinacrine 100 mg TID for 10 days.
        regimen.doses.QNC = local_doses(start_h, 8, 10, 100);


    %% --------------------------------------------------------------------
    %  Albendazole + nitroimidazole combinations
    %
    %  These are the important sensitivity schedules because the published
    %  refractory literature often aggregates "Alb + nitroimidazole" without
    %  cleanly separating MTZ+ABZ, TDZ+ABZ, exact doses, or durations.
    %  --------------------------------------------------------------------

    case 'MTZ5D_ABZ5D'
        regimen.doses.MTZ = local_doses(start_h, 8, 5, 400);
        regimen.doses.ABZ = local_doses(start_h, 24, 5, 400);

    case 'MTZ7D_ABZ7D'
        regimen.doses.MTZ = local_doses(start_h, 8, 7, 400);
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);

    case 'MTZ10D_ABZ7D'
        regimen.doses.MTZ = local_doses(start_h, 8, 10, 400);
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);

    case 'MTZ14D_ABZ7D'
        % Representative long-course combination used in current simulations.
        regimen.doses.MTZ = local_doses(start_h, 8, 14, 400);
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);

    case 'MTZ14D_ABZ14D'
        regimen.doses.MTZ = local_doses(start_h, 8, 14, 400);
        regimen.doses.ABZ = local_doses(start_h, 24, 14, 400);

    case 'TDZ2G_ABZ5D'
        regimen.doses.TDZ = local_doses(start_h, 24, 1, 2000);
        regimen.doses.TDZ = regimen.doses.TDZ(1);
        regimen.doses.ABZ = local_doses(start_h, 24, 5, 400);

    case 'TDZ2G_ABZ7D'
        regimen.doses.TDZ = local_doses(start_h, 24, 1, 2000);
        regimen.doses.TDZ = regimen.doses.TDZ(1);
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);

    case 'TDZ2G_ABZ14D'
        regimen.doses.TDZ = local_doses(start_h, 24, 1, 2000);
        regimen.doses.TDZ = regimen.doses.TDZ(1);
        regimen.doses.ABZ = local_doses(start_h, 24, 14, 400);


    %% --------------------------------------------------------------------
    %  Sequential treatment examples
    %
    %  These are not necessarily literature-matched. They are useful for
    %  model experiments asking whether a second class can rescue after a
    %  failed or partially suppressive first regimen.
    %  --------------------------------------------------------------------

    case 'MTZ5D_THEN_ABZ5D'
        regimen.doses.MTZ = local_doses(start_h, 8, 5, 400);
        regimen.doses.ABZ = local_doses(start_h + 24*5, 24, 5, 400);

    case 'MTZ7D_THEN_ABZ7D'
        regimen.doses.MTZ = local_doses(start_h, 8, 7, 400);
        regimen.doses.ABZ = local_doses(start_h + 24*7, 24, 7, 400);

    case 'ABZ5D_THEN_MTZ5D'
        regimen.doses.ABZ = local_doses(start_h, 24, 5, 400);
        regimen.doses.MTZ = local_doses(start_h + 24*5, 8, 5, 400);

    case 'ABZ7D_THEN_MTZ7D'
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);
        regimen.doses.MTZ = local_doses(start_h + 24*7, 8, 7, 400);

    case 'NTZ3D_THEN_QNC7D'
        regimen.doses.NTZ = local_doses(start_h, 12, 3, 500);
        regimen.doses.QNC = local_doses(start_h + 24*3, 8, 7, 100);

    case 'MTZ7D_ABZ7D_THEN_QNC7D'
        regimen.doses.MTZ = local_doses(start_h, 8, 7, 400);
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);
        regimen.doses.QNC = local_doses(start_h + 24*7, 8, 7, 100);

    case 'MTZ14D_ABZ7D_THEN_QNC7D'
        regimen.doses.MTZ = local_doses(start_h, 8, 14, 400);
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);
        regimen.doses.QNC = local_doses(start_h + 24*14, 8, 7, 100);


    %% --------------------------------------------------------------------
    %  Rescue / exploratory combinations
    %  --------------------------------------------------------------------

    case 'ABZ7D_QNC7D'
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);
        regimen.doses.QNC = local_doses(start_h, 8, 7, 100);

    case 'NTZ3D_ABZ7D'
        regimen.doses.NTZ = local_doses(start_h, 12, 3, 500);
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);

    case 'NTZ7D_ABZ7D'
        regimen.doses.NTZ = local_doses(start_h, 12, 7, 500);
        regimen.doses.ABZ = local_doses(start_h, 24, 7, 400);

    case 'NTZ3D_QNC7D'
        regimen.doses.NTZ = local_doses(start_h, 12, 3, 500);
        regimen.doses.QNC = local_doses(start_h, 8, 7, 100);

    otherwise
        error('Unknown regimen "%s". Normalised key was "%s".', which, key);
end

end

% -------------------------------------------------------------------------
function D = local_doses(start_h, interval_h, duration_days, amount_mg)
% Build repeated dose structs.
%
% start_h       first dose time
% interval_h    time between doses
% duration_days dosing duration
% amount_mg     dose amount

stop_h = start_h + 24*duration_days;
t = start_h:interval_h:(stop_h - 1e-9);

D = arrayfun(@(ti) struct('time_h', ti, 'amount_mg', amount_mg), t);
end

%{
% Regimen names: 

MTZ5D_ABZ5D
MTZ7D_ABZ7D
MTZ10D_ABZ7D
MTZ14D_ABZ7D
MTZ14D_ABZ14D

TDZ2G_ABZ5D
TDZ2G_ABZ7D
TDZ2G_ABZ14D

MTZ5D_THEN_ABZ5D
MTZ7D_THEN_ABZ7D
ABZ5D_THEN_MTZ5D
ABZ7D_THEN_MTZ7D

NTZ3D
NTZ7D
QNC5D
QNC7D
QNC10D

NTZ / ABZ / QNC exploratory combinations
%}