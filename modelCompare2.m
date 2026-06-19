function summary_table = modelCompare2_updated(t_end_days, make_plots)
% MODELCOMPARE2_UPDATED
% Deterministic Giardia protocol comparison with NTZ and QNC included.
%
% This version is designed for pre-Monte-Carlo checking: it runs the model
% with the current baseline parameters, plus a small set of local literature-
% informed refinements for deterministic comparison.
%
% Key features:
%   - includes untreated, MTZ5D, ABZ5D, NTZ3D, TDZ2G_ONCE,
%     MTZ14D_ABZ7D, and QNC7D;
%   - compares baseline vs a nitroimidazole-reduced scenario;
%   - reduced susceptibility is applied to MTZ and TDZ only here;
%   - ABZ, NTZ, and QNC are left as separate PD blocks;
%   - outputs clearance, rebound, durable clearance, min/day30/final burden;
%   - saves a CSV summary table.
%
% Usage:
%   summary_table = modelCompare2_updated();
%   summary_table = modelCompare2_updated(30, true);
%
% Requirements:
%   giardia_params.m
%   giardia_regimens.m
%   giardia_simulate.m
%
% Note:
%   If you want this to replace the old file, rename it to modelCompare2.m.

if nargin < 1 || isempty(t_end_days)
    t_end_days = 30;
end
if nargin < 2 || isempty(make_plots)
    make_plots = true;
end

p0 = giardia_params();
p0 = local_refine_parameters_for_deterministic_check(p0);

tvec = 0:0.1:(24*t_end_days);  % hours; 0.1 h gives smooth deterministic curves

% Focused deterministic comparison, including nitazoxanide and quinacrine.
regimen_keys = { ...
    'NONE', ...
    'MTZ5D', ...
    'ABZ5D', ...
    'NTZ3D', ...
    'TDZ2G_ONCE', ...
    'MTZ14D_ABZ7D', ...
    'QNC7D'};

regimen_labels = { ...
    'Untreated', ...
    'MTZ5D', ...
    'ABZ5D', ...
    'NTZ3D', ...
    'TDZ2G', ...
    'MTZ14D+ABZ7D', ...
    'QNC7D'};

% Scenario names. The reduced-effect scenario is deliberately nitroimidazole-
% specific; it should not automatically penalise ABZ, NTZ, or QNC.
scenario_labels = {'baseline', 'nitroimidazole-reduced'};
scenario_nitro_reduced = [false, true];

n_rows = numel(regimen_keys) * numel(scenario_labels);

Scenario = cell(n_rows,1);
Protocol = cell(n_rows,1);
Cleared = false(n_rows,1);
Clearance_time_days = NaN(n_rows,1);
Rebounded = false(n_rows,1);
Durable_clearance = false(n_rows,1);
Min_T = NaN(n_rows,1);
Day30_T = NaN(n_rows,1);
Final_T = NaN(n_rows,1);

outs = cell(numel(scenario_labels), numel(regimen_keys));
row = 0;

for s = 1:numel(scenario_labels)
    p = p0;

    if scenario_nitro_reduced(s)
        p = local_apply_nitroimidazole_reduced_effect(p);
    end

    for r = 1:numel(regimen_keys)
        row = row + 1;

        reg = local_make_regimen(regimen_keys{r});

        % opts.resistant is kept false because this script applies the
        % nitroimidazole-reduced scenario directly to MTZ/TDZ parameters.
        % This avoids accidentally reducing ABZ, NTZ, or QNC in older
        % versions of giardia_simulate.m.
        opts = struct('useMicrobiome', true, 'resistant', false);
        out = giardia_simulate(tvec, reg, p, opts);
        outs{s,r} = out;

        Scenario{row} = scenario_labels{s};
        Protocol{row} = regimen_labels{r};

        if ~isnan(out.clearance_time_h)
            Clearance_time_days(row) = out.clearance_time_h / 24;
            Cleared(row) = true;
        end

        Rebounded(row) = logical(out.rebounded);
        Durable_clearance(row) = Cleared(row) && ~Rebounded(row);
        Min_T(row) = min(out.T);
        Final_T(row) = out.T(end);

        if t_end_days >= 30
            Day30_T(row) = interp1(out.tvec(:)/24, out.T(:), 30, 'linear');
        end
    end
end

summary_table = table(Scenario, Protocol, Cleared, Clearance_time_days, ...
    Rebounded, Durable_clearance, Min_T, Day30_T, Final_T);

disp(' ')
disp('=== Deterministic protocol summary: baseline and nitroimidazole-reduced scenarios ===')
disp(summary_table)

try
    writetable(summary_table, 'modelCompare2_deterministic_summary.csv');
    fprintf('\nSaved: modelCompare2_deterministic_summary.csv\n');
catch ME
    warning('Could not write CSV file: %s', ME.message);
end

if make_plots
    local_plot_scenarios(outs, regimen_labels, scenario_labels, p0, t_end_days);
end

end

% -------------------------------------------------------------------------
function p = local_refine_parameters_for_deterministic_check(p)
% Local deterministic refinements.
% These are intentionally conservative and transparent. The parameter file
% remains the canonical source; this helper just makes the comparison script
% less dependent on older placeholders.

% Tinidazole: mean half-life after 2 g oral dosing is around 13 h.
p.TDZ.halflife_h = 13.2;
p.TDZ.kel = log(2) / p.TDZ.halflife_h;

% Nitazoxanide/tizoxanide proxy. Keep close to the current parameter file:
% short systemic half-life, with the main uncertainty living in local gut
% exposure and PD rather than systemic persistence.
p.NTZ.halflife_h = 1.3;
p.NTZ.kel = log(2) / p.NTZ.halflife_h;

% Quinacrine is used here as a legacy/rescue drug with long persistence.
% This is still exploratory, but 24 h was too short for a deterministic
% rescue-comparison placeholder.
p.QNC.halflife_h = 120.0;
p.QNC.kel = log(2) / p.QNC.halflife_h;

% Keep the plasma contribution explicit if present in the parameter file.
% giardia_pk.m may still hard-code this unless separately patched.
if ~isfield(p, 'exposure')
    p.exposure.alpha_plasma = 0.2;
elseif ~isfield(p.exposure, 'alpha_plasma')
    p.exposure.alpha_plasma = 0.2;
end

end

% -------------------------------------------------------------------------
function p = local_apply_nitroimidazole_reduced_effect(p)
% Apply reduced susceptibility to nitroimidazoles only.
% This mirrors the current Strategy B interpretation:
% refractory phenotype mainly reduces MTZ/TDZ effect, while ABZ, NTZ, and
% QNC remain separate PD blocks.

scale = p.resistance.scale;

p.MTZ.Emax = scale * p.MTZ.Emax;
p.TDZ.Emax = scale * p.TDZ.Emax;

% Optional concentration-scale penalty. Leave commented for now because the
% previous model convention used Emax scaling only.
% p.MTZ.MIC = p.MTZ.MIC / scale;
% p.TDZ.MIC = p.TDZ.MIC / scale;

end

% -------------------------------------------------------------------------
function reg = local_make_regimen(key)
% Create regimen, including an empty no-treatment regimen.
if strcmpi(key, 'NONE')
    reg.name = 'NONE';
    reg.doses = struct();
else
    reg = giardia_regimens(key, 0);
end
end

% -------------------------------------------------------------------------
function local_plot_scenarios(outs, regimen_labels, scenario_labels, p, t_end_days)

figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

floor_T = 1e-15;

for s = 1:numel(scenario_labels)
    nexttile;
    hold on;

    h = gobjects(numel(regimen_labels),1);

    for r = 1:numel(regimen_labels)
        out = outs{s,r};
        y = max(out.T(:), floor_T);
        h(r) = plot(out.tvec(:)/24, y, 'LineWidth', 2);
    end

    set(gca, 'YScale', 'log');
    yline(p.clearance_threshold, 'k--', 'Clearance threshold');

    xlabel('Time (days)');
    ylabel('Trophozoites (billions, log scale)');
    title(sprintf('%s scenario', scenario_labels{s}));
    xlim([0 t_end_days]);
    ylim([floor_T 20]);

    legend(h, regimen_labels, 'Location', 'eastoutside');
    box on;
    grid on;
    grid minor;
end

end
