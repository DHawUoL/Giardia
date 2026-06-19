function p = giardia_params()
% GIARDIA_PARAMS  Parameter set, uncertainty metadata, and calibration targets.
%
% This is a literature-audit draft. The fields used directly by the simulator
% are still simple numeric values, e.g. p.rT, p.MTZ.Emax, p.ABZ.F_duod.
% Additional metadata are stored in:
%   p.refs        source handles used in notes/tables
%   p.blocks      block names for structured uncertainty analysis
%   p.u           parameter uncertainty/status metadata
%   p.calibration clinical calibration/validation targets
%   p.populations susceptible/refractory population modifiers
%
% Units:
%   time: hours
%   T, C: 1e9 organisms/cysts unless otherwise stated
%   concentrations: ug/L in the current PK/PD code
%
% Model states:
%   T = trophozoites in the upper small intestine
%   C = cyst burden/output proxy in the modelled compartment
%   M = protective microbiome proxy, scaled to [0,1]
%
% Current model equations, treatment-free except for kdrug(t):
%   dT/dt = rT*T*(1 - T/K) - k_encyst*T - k_I*T ...
%           - k_micro*M*T - kdrug(t)*T
%   dC/dt = k_encyst*T - k_c_decay*C
%   dM/dt = rM*M*(1 - M) - dys*M
%
% Current PK structure in giardia_pk.m:
%   oral dose -> stomach -> duodenum -> plasma
%   Ceff = Cduod + alpha*Cplasma
%
% IMPORTANT:
%   Many biological and PD parameters are effective/calibrated quantities.
%   They should not be described as direct literature measurements unless
%   p.u.<parameter>.status says so.

%% ------------------------------------------------------------------------
%  Source handles
%  Full references should live in the manuscript/spreadsheet. These handles
%  are for keeping the MATLAB metadata readable.
%% ------------------------------------------------------------------------

p.refs.FinkSinger2020 = 'Fink & Singer 2020: in-vitro Giardia growth/doubling-time context';
p.refs.Adam2021       = 'Adam 2021: Giardia biology/pathogenesis review; drug half-life table/context';
p.refs.Thomas2021     = 'Thomas et al. 2021: encystation staging/variability context';
p.refs.Karabay2004    = 'Karabay et al. 2004: adult ABZ5D vs MTZ5D clinical comparison';
p.refs.Morch2008      = 'Morch/Kristine et al. 2008: refractory giardiasis treatment ladder';
p.refs.Mattila1983    = 'Mattila et al. 1983: comparative MTZ/TDZ pharmacokinetics';
p.refs.Labels         = 'Drug labels / prescribing information';
p.refs.Block5Reviews  = 'Refractory giardiasis reviews: context, regimens, scenario justification';
p.refs.Placeholder    = 'Placeholder pending literature check or calibration';
p.refs.Calibrated     = 'Calibrated/effective model parameter, not direct literature estimate';

%% ------------------------------------------------------------------------
%  Structured uncertainty blocks
%% ------------------------------------------------------------------------

p.blocks.parasite_natural_history = {'rT','K','k_encyst','k_c_decay','T0','C0'};
p.blocks.host_microbiome          = {'k_I','k_micro','M0','rM','dys'};
p.blocks.shared_gut_physiology    = {'k_GE','k_tr_duod','V_duod_L'};
p.blocks.systemic_PK              = {'halflife_h','kel','Vd_L'};
p.blocks.gut_site_exposure        = {'F_duod','V_duod_L','alpha_plasma'};
p.blocks.PD_nitroimidazole        = {'MTZ.Emax','MTZ.MIC','MTZ.h','TDZ.Emax','TDZ.MIC','TDZ.h'};
p.blocks.PD_benzimidazole         = {'ABZ.Emax','ABZ.MIC','ABZ.h'};
p.blocks.PD_alternative           = {'NTZ.Emax','NTZ.MIC','NTZ.h','QNC.Emax','QNC.MIC','QNC.h'};
p.blocks.scenario                 = {'resistance.frac','resistance.scale','adherence.p_miss'};

%% ------------------------------------------------------------------------
%  Growth / interaction parameters
%% ------------------------------------------------------------------------

p.rT       = 0.04;      % effective trophozoite growth rate, per h
p.K        = 20;        % effective carrying capacity, billions
p.k_encyst = 0.002;     % encystation / cyst-production rate, per h
p.k_I      = 0.001;     % baseline host clearance, per h
p.k_micro  = 0.02;      % microbiome suppression coefficient
p.k_c_decay= 0.02;      % cyst loss/removal rate, per h
p.dys      = 0.005;     % dysbiosis dampening of M, per h
p.rM       = 0.02;      % microbiome recovery rate, per h

%% Initial conditions
p.T0       = 0.5;       % initial trophozoites, billions
p.C0       = 0.1;       % initial cyst burden/output proxy, billions
p.M0       = 0.3;       % initial microbiome protection, [0,1]

%% Operational clearance definition
p.clearance_threshold = 1e-3;   % billions (~1e6 trophozoites)
p.clearance_hold_h    = 72;     % below threshold for 72 h

%% ------------------------------------------------------------------------
%  Systemic PK parameters
%  These are the best candidates to fix in the primary analysis and vary
%  only narrowly. kel is always derived from half-life.
%% ------------------------------------------------------------------------

% Metronidazole (MTZ)
p.MTZ.halflife_h = 8.0;
p.MTZ.kel        = log(2)/p.MTZ.halflife_h;
p.MTZ.Vd_L       = 50;
p.MTZ.f_lumen    = NaN; % deprecated older-model parameter; retained only for compatibility

% Tinidazole (TDZ)
p.TDZ.halflife_h = 13.0;
p.TDZ.kel        = log(2)/p.TDZ.halflife_h;
p.TDZ.Vd_L       = 50;
p.TDZ.f_lumen    = NaN; % deprecated

% Albendazole (ABZ) -> active metabolite sulfoxide proxy
p.ABZ.halflife_h = 10.0;
p.ABZ.kel        = log(2)/p.ABZ.halflife_h;
p.ABZ.Vd_L       = 20;
p.ABZ.f_lumen    = NaN; % deprecated

% Nitazoxanide (NTZ) -> tizoxanide-effective exposure proxy
% Placeholder until NTZ/tizoxanide PK/PD are checked carefully.
p.NTZ.halflife_h = 1.3;
p.NTZ.kel        = log(2)/p.NTZ.halflife_h;
p.NTZ.Vd_L       = 20;
p.NTZ.f_lumen    = NaN; % deprecated

% Quinacrine (QNC) legacy/refractory rescue proxy
p.QNC.halflife_h = 24.0;
p.QNC.kel        = log(2)/p.QNC.halflife_h;
p.QNC.Vd_L       = 150;
p.QNC.f_lumen    = NaN; % deprecated

%% ------------------------------------------------------------------------
%  Pharmacodynamic kill parameters
%  These remain calibrated/effective unless later fitted to susceptibility
%  or clinical response data.
%% ------------------------------------------------------------------------

p.MTZ.Emax = 0.08; p.MTZ.MIC = 20; p.MTZ.h = 1.2;  % tuned refractory-prone baseline
p.TDZ.Emax = 0.35; p.TDZ.MIC = 4;  p.TDZ.h = 1.2;
p.ABZ.Emax = 0.15; p.ABZ.MIC = 1;  p.ABZ.h = 1.5;
p.NTZ.Emax = 0.20; p.NTZ.MIC = 5;  p.NTZ.h = 1.2;  % placeholder
p.QNC.Emax = 0.25; p.QNC.MIC = 3;  p.QNC.h = 1.2;

%% Adherence and resistance scenario parameters
p.adherence.p_miss = 0.05;
p.resistance.frac  = 0.15;
p.resistance.scale = 0.3;  % Emax -> scale*Emax when resistant=true

%% ------------------------------------------------------------------------
%  Minimal stomach -> duodenum -> plasma PK parameters
%% ------------------------------------------------------------------------

% Global plasma contribution to effective exposure.
% NOTE: giardia_pk.m currently hard-codes alpha=0.2 unless updated to read this.
p.exposure.alpha_plasma = 0.2;

% Metronidazole
p.MTZ.k_GE       = 1.2;
p.MTZ.k_abs_duod = 1.5;
p.MTZ.k_tr_duod  = 0.8;
p.MTZ.V_duod_L   = 0.25;
p.MTZ.F_duod     = 0.01;

% Tinidazole
p.TDZ.k_GE       = 1.0;
p.TDZ.k_abs_duod = 1.2;
p.TDZ.k_tr_duod  = 0.8;
p.TDZ.V_duod_L   = 0.25;
p.TDZ.F_duod     = 0.01;

% Albendazole
p.ABZ.k_GE       = 1.0;
p.ABZ.k_abs_duod = 0.3;
p.ABZ.k_tr_duod  = 0.8;
p.ABZ.V_duod_L   = 0.25;
p.ABZ.F_duod     = 0.003;

% Nitazoxanide / tizoxanide proxy
p.NTZ.k_GE       = 1.0;
p.NTZ.k_abs_duod = 0.8;
p.NTZ.k_tr_duod  = 0.8;
p.NTZ.V_duod_L   = 0.25;
p.NTZ.F_duod     = 0.01;

% Quinacrine
p.QNC.k_GE       = 1.0;
p.QNC.k_abs_duod = 0.6;
p.QNC.k_tr_duod  = 0.8;
p.QNC.V_duod_L   = 0.25;
p.QNC.F_duod     = 0.01;

%% ------------------------------------------------------------------------
%  Uncertainty metadata
%  range convention:
%    triangular: [min mode max]
%    lognormal_cv: [CV] or [lo mode hi] if a bounded display range is useful
%    fixed: []
%  These are draft ranges. status distinguishes literature-informed values
%  from calibrated/phenomenological/placeholders.
%% ------------------------------------------------------------------------

p.u.rT        = U('parasite_natural_history','literature_informed_effective','triangular',[0.025 0.04 0.09],'FinkSinger2020','In-vitro doubling times inform scale; baseline is slower effective in-host growth.');
p.u.K         = U('parasite_natural_history','calibrated_effective','lognormal_cv',[0.5],'Calibrated','Effective carrying capacity; calibrated to untreated persistence.');
p.u.k_encyst  = U('parasite_natural_history','weakly_informed_phenomenological','lognormal_cv',[0.75],'Thomas2021','Encystation biology known, but no clean in-vivo rate.');
p.u.k_c_decay = U('parasite_natural_history','phenomenological','lognormal_cv',[0.75],'Placeholder','Cyst loss/output removal parameter, not directly estimated.');
p.u.T0        = U('parasite_natural_history','initial_condition','lognormal_cv',[0.75],'Placeholder','Initial trophozoite burden; vary as initial-condition sensitivity.');
p.u.C0        = U('parasite_natural_history','initial_condition','lognormal_cv',[0.75],'Placeholder','Initial cyst burden/output proxy; correlated with T0.');

p.u.k_I       = U('host_microbiome','phenomenological','lognormal_cv',[0.75],'Placeholder','Baseline host clearance; not directly estimated.');
p.u.k_micro   = U('host_microbiome','phenomenological','lognormal_cv',[0.75],'Placeholder','Microbiome suppression coefficient; scenario/sensitivity parameter.');
p.u.M0        = U('host_microbiome','scenario_bounded','bounded_normal',[0 0.3 1],'Placeholder','Initial microbiome protection.');
p.u.rM        = U('host_microbiome','phenomenological','lognormal_cv',[0.75],'Placeholder','Microbiome recovery rate.');
p.u.dys       = U('host_microbiome','phenomenological','lognormal_cv',[0.75],'Placeholder','Dysbiosis dampening parameter.');

p.u.clearance_threshold = U('observation','definition','fixed',[],'Calibrated','Operational reporting threshold; sensitivity possible.');
p.u.clearance_hold_h    = U('observation','definition','fixed',[],'Calibrated','Operational clearance hold time; sensitivity possible.');

% Systemic PK: fixed/narrow primary analysis
p.u.MTZ.halflife_h = U('systemic_PK','fixed_from_label_or_PK','triangular',[7.5 8.0 9.0],'Mattila1983/Labels','MTZ half-life is relatively well constrained.');
p.u.TDZ.halflife_h = U('systemic_PK','fixed_from_label_or_PK','triangular',[12.0 13.5 15.0],'Mattila1983/Labels','TDZ half-life relatively well constrained.');
p.u.ABZ.halflife_h = U('systemic_PK','fixed_from_label_or_PK','triangular',[8.0 10.0 12.0],'Labels','Active albendazole sulfoxide half-life proxy.');
p.u.NTZ.halflife_h = U('systemic_PK','placeholder_needs_check','triangular',[1.0 1.3 2.0],'Labels','NTZ/tizoxanide value needs confirmation before use.');
p.u.QNC.halflife_h = U('systemic_PK','placeholder_needs_check','triangular',[12.0 24.0 48.0],'Placeholder','QNC systemic PK requires further checking.');

for nm = {'MTZ','TDZ','ABZ','NTZ','QNC'}
    drug = nm{1};
    p.u.(drug).kel = U('systemic_PK','derived','fixed',[],'derived','kel = log(2)/halflife_h; do not sample independently.');
    p.u.(drug).Vd_L = U('systemic_PK','approximate','lognormal_cv',[0.25],'Labels','Fix or vary narrowly/moderately; not main uncertainty bottleneck.');
    p.u.(drug).k_GE = U('shared_gut_physiology','effective_uncertain','lognormal_cv',[0.3],'Placeholder','Shared gastric-emptying/gut physiology uncertainty.');
    p.u.(drug).k_tr_duod = U('shared_gut_physiology','effective_uncertain','lognormal_cv',[0.3],'Placeholder','Shared duodenal transit uncertainty.');
    p.u.(drug).V_duod_L = U('gut_site_exposure','effective_calibrated','lognormal_cv',[0.5],'Calibrated','Effective duodenal volume, sensitivity parameter.');
    p.u.(drug).k_abs_duod = U('drug_specific_PK','effective_uncertain','lognormal_cv',[0.5],'Calibrated','Effective duodenal absorption into plasma.');
    p.u.(drug).F_duod = U('gut_site_exposure','effective_calibrated','lognormal_cv',[0.75],'Calibrated','Effective local availability; key sensitivity/calibration parameter.');
end

p.u.exposure.alpha_plasma = U('gut_site_exposure','effective_calibrated','triangular',[0.0 0.2 0.5],'Calibrated','Weight of plasma contribution to effective exposure.');

% PD metadata
p.u.MTZ.Emax = U('PD_nitroimidazole','calibrated_effective','lognormal_cv',[0.75],'Karabay2004/Morch2008','Tune/calibrate via susceptible and refractory clinical targets.');
p.u.MTZ.MIC  = U('PD_nitroimidazole','calibrated_effective','lognormal_cv',[0.75],'susceptibility_lit_needed','In-vivo concentration scale; not direct literature MIC yet.');
p.u.MTZ.h    = U('PD_nitroimidazole','weakly_constrained','triangular',[0.8 1.2 2.0],'Placeholder','Hill coefficient sensitivity.');

p.u.TDZ.Emax = U('PD_nitroimidazole','calibrated_effective','lognormal_cv',[0.75],'Morch2008/susceptibility_lit_needed','TDZ linked to nitroimidazole susceptibility block.');
p.u.TDZ.MIC  = U('PD_nitroimidazole','calibrated_effective','lognormal_cv',[0.75],'susceptibility_lit_needed','TDZ concentration scale.');
p.u.TDZ.h    = U('PD_nitroimidazole','weakly_constrained','triangular',[0.8 1.2 2.0],'Placeholder','Hill coefficient sensitivity.');

p.u.ABZ.Emax = U('PD_benzimidazole','calibrated_effective','lognormal_cv',[0.5],'Karabay2004/Morch2008','ABZ constrained by susceptible comparison and refractory combination therapy.');
p.u.ABZ.MIC  = U('PD_benzimidazole','calibrated_effective','lognormal_cv',[0.5],'susceptibility_lit_needed','ABZ concentration scale.');
p.u.ABZ.h    = U('PD_benzimidazole','weakly_constrained','triangular',[0.8 1.5 2.5],'Placeholder','Hill coefficient sensitivity.');

p.u.NTZ.Emax = U('PD_alternative','placeholder_needs_check','lognormal_cv',[0.75],'Labels/susceptibility_lit_needed','NTZ monotherapy to add; values are placeholders.');
p.u.NTZ.MIC  = U('PD_alternative','placeholder_needs_check','lognormal_cv',[0.75],'susceptibility_lit_needed','NTZ concentration scale placeholder.');
p.u.NTZ.h    = U('PD_alternative','placeholder_needs_check','triangular',[0.8 1.2 2.0],'Placeholder','Hill coefficient placeholder.');

p.u.QNC.Emax = U('PD_alternative','calibrated_effective_or_placeholder','lognormal_cv',[0.75],'Morch2008/QNC_lit_needed','QNC rescue efficacy; small-denominator calibration target.');
p.u.QNC.MIC  = U('PD_alternative','placeholder_needs_check','lognormal_cv',[0.75],'susceptibility_lit_needed','QNC concentration scale placeholder.');
p.u.QNC.h    = U('PD_alternative','weakly_constrained','triangular',[0.8 1.2 2.0],'Placeholder','Hill coefficient sensitivity.');

p.u.adherence.p_miss = U('scenario','scenario','fixed_or_scenario',[],'Block5Reviews','Behavioural/protocol scenario rather than biological covariance.');
p.u.resistance.frac  = U('scenario','scenario','fixed_or_scenario',[],'Morch2008','Defines mixture/refractory fraction; not primary biological uncertainty.');
p.u.resistance.scale = U('scenario','scenario','fixed_or_scenario',[],'Morch2008/susceptibility_lit_needed','Reduced-susceptibility efficacy multiplier.');

%% ------------------------------------------------------------------------
%  Clinical calibration / validation targets
%% ------------------------------------------------------------------------

% Susceptible/non-refractory adult response: Karabay et al.
p.calibration.susceptible.Karabay2004.note = 'Adult ordinary giardiasis; use for susceptible/non-refractory validation/calibration, not refractory scenario.';
p.calibration.susceptible.Karabay2004.MTZ5D.n = 29;
p.calibration.susceptible.Karabay2004.MTZ5D.day7_negative = 29;
p.calibration.susceptible.Karabay2004.MTZ5D.day15_negative = 29;
p.calibration.susceptible.Karabay2004.ABZ5D.n = 28;
p.calibration.susceptible.Karabay2004.ABZ5D.day7_negative = 28;
p.calibration.susceptible.Karabay2004.ABZ5D.day15_negative = 27;
p.calibration.susceptible.Karabay2004.endpoint = 'faecal examination negative at follow-up';

% Refractory treatment ladder: Morch/Kristine et al.
p.calibration.refractory.Morch2008.note = 'Metronidazole-refractory giardiasis treatment ladder after Bergen outbreak.';
p.calibration.refractory.Morch2008.MTZ_ABZ.n = 38;
p.calibration.refractory.Morch2008.MTZ_ABZ.effective = 30;
p.calibration.refractory.Morch2008.MTZ_ABZ.failed = 8;
p.calibration.refractory.Morch2008.PAROMOMYCIN.n = 6;
p.calibration.refractory.Morch2008.PAROMOMYCIN.effective = 3;
p.calibration.refractory.Morch2008.QNC_MTZ.n = 3;
p.calibration.refractory.Morch2008.QNC_MTZ.effective = 3;
p.calibration.refractory.Morch2008.endpoint = 'observed effective response; extract exact retest timing/method from tables when available';
p.calibration.refractory.Morch2008.to_extract = {'days_post_treatment_to_retest','diagnostic_method','positivity_definition','symptom_vs_parasitological_cure'};

%% ------------------------------------------------------------------------
%  Population modifiers: Strategy B
%% ------------------------------------------------------------------------

p.populations.susceptible.description = 'Ordinary/non-refractory giardiasis. Calibrate/validate with Karabay-type MTZ/ABZ high response.';
p.populations.susceptible.nitro_Emax_multiplier = 1.0;
p.populations.susceptible.nitro_MIC_multiplier  = 1.0;
p.populations.susceptible.ABZ_Emax_multiplier   = 1.0;
p.populations.susceptible.ABZ_MIC_multiplier    = 1.0;

p.populations.refractory.description = 'Reduced-susceptibility/refractory phenotype. Calibrate with Morch treatment ladder.';
p.populations.refractory.nitro_Emax_multiplier = p.resistance.scale;
p.populations.refractory.nitro_MIC_multiplier  = 1/p.resistance.scale;
p.populations.refractory.ABZ_Emax_multiplier   = 1.0;  % separate block; do not assume global resistance
p.populations.refractory.ABZ_MIC_multiplier    = 1.0;
p.populations.refractory.note = 'Draft: refractory primarily shifts nitroimidazole susceptibility; ABZ remains partly independent.';

end

% -------------------------------------------------------------------------
function out = U(block, status, dist, range, ref, note)
% U  Small helper for parameter uncertainty metadata.
out = struct();
out.block  = block;
out.status = status;
out.dist   = dist;
out.range  = range;
out.ref    = ref;
out.note   = note;
end
