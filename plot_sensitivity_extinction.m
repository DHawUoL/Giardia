C = lines(7);

qnc_file = 'mcmc_outputs/pooled_QNC7D_4chains_currentSimulator.mat';
mtz_file = 'mcmc_outputs/pooled_MTZ14D_ABZ7D_4chains_currentSimulator.mat';

cases = struct([]);

% -------------------------------------------------------------------------
% QNC
% -------------------------------------------------------------------------
cases(1).chain_file = qnc_file;
cases(1).regimen = 'QNC7D';
cases(1).scenario = 'refractory';
cases(1).label = 'QNC7D, residual persistence';
cases(1).absorbExtinction = false;
cases(1).extinction_threshold = NaN;
cases(1).line_style = '-';
cases(1).color = C(1,:);

cases(2).chain_file = qnc_file;
cases(2).regimen = 'QNC7D';
cases(2).scenario = 'refractory';
cases(2).label = 'QNC7D, AE, T_{ext}=10^{-9}';
cases(2).absorbExtinction = true;
cases(2).extinction_threshold = 1e-9;
cases(2).line_style = '--';
cases(2).color = C(1,:);

cases(3).chain_file = qnc_file;
cases(3).regimen = 'QNC7D';
cases(3).scenario = 'refractory';
cases(3).label = 'QNC7D, AE, T_{ext}=10^{-12}';
cases(3).absorbExtinction = true;
cases(3).extinction_threshold = 1e-12;
cases(3).line_style = ':';
cases(3).color = C(1,:);

% -------------------------------------------------------------------------
% MTZ14D + ABZ7D
% -------------------------------------------------------------------------
cases(4).chain_file = mtz_file;
cases(4).regimen = 'MTZ14D_ABZ7D';
cases(4).scenario = 'refractory';
cases(4).label = 'MTZ14D+ABZ7D, residual persistence';
cases(4).absorbExtinction = false;
cases(4).extinction_threshold = NaN;
cases(4).line_style = '-';
cases(4).color = C(7,:);

cases(5).chain_file = mtz_file;
cases(5).regimen = 'MTZ14D_ABZ7D';
cases(5).scenario = 'refractory';
cases(5).label = 'MTZ14D+ABZ7D, AE, T_{ext}=10^{-9}';
cases(5).absorbExtinction = true;
cases(5).extinction_threshold = 1e-9;
cases(5).line_style = '--';
cases(5).color = C(7,:);

cases(6).chain_file = mtz_file;
cases(6).regimen = 'MTZ14D_ABZ7D';
cases(6).scenario = 'refractory';
cases(6).label = 'MTZ14D+ABZ7D, AE, T_{ext}=10^{-12}';
cases(6).absorbExtinction = true;
cases(6).extinction_threshold = 1e-12;
cases(6).line_style = ':';
cases(6).color = C(7,:);

summary = plot_outcome_probability_curves_cases(cases);