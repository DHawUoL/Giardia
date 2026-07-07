function chain_pool = pool_chains_from_files(chain_files, burn_frac, out_file)
%POOL_CHAINS_FROM_FILES Pool MCMC chains loaded from .mat files.
%
% Usage:
%   chain_pool = pool_chains_from_files(chain_files, 0.5);
%   chain_pool = pool_chains_from_files(chain_files, 0.5, out_file);
%
% Inputs:
%   chain_files : cell array / string array of .mat filenames, OR a wildcard pattern
%                 e.g. 'mcmc_outputs/giardia_mcmc_qnc_refractory_50_54_mcmc*_vp1001_chain.mat'
%
%   burn_frac   : fraction of each chain to discard before pooling.
%                 Default 0.5.
%
%   out_file    : optional .mat file to save pooled chain into.
%
% Output:
%   chain_pool  : pooled chain struct.

if nargin < 2 || isempty(burn_frac)
    burn_frac = 0.5;
end

if nargin < 3
    out_file = '';
end

% -------------------------------------------------------------------------
% Resolve file list
% -------------------------------------------------------------------------
if ischar(chain_files) || isstring(chain_files)
    chain_files = string(chain_files);

    % If a wildcard/pattern was supplied, expand it.
    if contains(chain_files, "*")
        D = dir(chain_files);
        if isempty(D)
            error('No files matched pattern: %s', chain_files);
        end

        folder = string({D.folder})';
        names  = string({D.name})';
        chain_files = fullfile(folder, names);

    else
        chain_files = cellstr(chain_files);
    end
end

chain_files = cellstr(chain_files);
n_chains = numel(chain_files);

if n_chains < 1
    error('No chain files supplied.');
end

% -------------------------------------------------------------------------
% Load chains
% -------------------------------------------------------------------------
chains = cell(n_chains,1);

for c = 1:n_chains
    this_file = chain_files{c};

    if ~isfile(this_file)
        error('File not found: %s', this_file);
    end

    S = load(this_file);

    if isfield(S, 'chain')
        ch = S.chain;
    else
        % Fall back: find the first struct with a theta field.
        fn = fieldnames(S);
        found = false;

        for k = 1:numel(fn)
            candidate = S.(fn{k});
            if isstruct(candidate) && isfield(candidate, 'theta')
                ch = candidate;
                found = true;
                break
            end
        end

        if ~found
            error('No chain struct found in file: %s', this_file);
        end
    end

    required_fields = {'theta','logpost','accepted','param_names'};
    for r = 1:numel(required_fields)
        if ~isfield(ch, required_fields{r})
            error('File %s contains a chain but is missing field: %s', ...
                this_file, required_fields{r});
        end
    end

    chains{c} = ch;
end

% -------------------------------------------------------------------------
% Check compatibility
% -------------------------------------------------------------------------
param_names_1 = chains{1}.param_names;

for c = 2:n_chains
    if numel(chains{c}.param_names) ~= numel(param_names_1) || ...
            any(~strcmp(chains{c}.param_names(:), param_names_1(:)))

        fprintf('\nChain 1 parameters:\n')
        disp(param_names_1(:))

        fprintf('\nChain %d parameters:\n', c)
        disp(chains{c}.param_names(:))

        error('Parameter names differ across chains. Do not pool these as one posterior.');
    end
end

% Optional warning if scenario names differ
if isfield(chains{1}, 'scenario_name')
    scenario_1 = string(chains{1}.scenario_name);

    for c = 2:n_chains
        if isfield(chains{c}, 'scenario_name')
            if string(chains{c}.scenario_name) ~= scenario_1
                warning('Scenario names differ across chains. Check this is intentional.');
            end
        end
    end
end

% -------------------------------------------------------------------------
% Pool post-burn samples
% -------------------------------------------------------------------------
chain_pool = chains{1};

Theta_pool = [];
LogPost_pool = [];
Accepted_pool = [];

chain_source = table( ...
    strings(n_chains,1), ...
    zeros(n_chains,1), ...
    zeros(n_chains,1), ...
    zeros(n_chains,1), ...
    'VariableNames', {'file','n_total','burn_index','n_kept'} ...
);

for c = 1:n_chains
    ch = chains{c};

    n_total = size(ch.theta,1);
    burn = floor(burn_frac * n_total) + 1;

    if burn > n_total
        error('Burn-in leaves no samples for chain %d.', c);
    end

    Theta_pool    = [Theta_pool; ch.theta(burn:end,:)];
    LogPost_pool  = [LogPost_pool; ch.logpost(burn:end)];
    Accepted_pool = [Accepted_pool; ch.accepted(burn:end)];

    chain_source.file(c) = string(chain_files{c});
    chain_source.n_total(c) = n_total;
    chain_source.burn_index(c) = burn;
    chain_source.n_kept(c) = n_total - burn + 1;
end

chain_pool.theta = Theta_pool;
chain_pool.logpost = LogPost_pool;
chain_pool.accepted = Accepted_pool;
chain_pool.accept_rate = mean(Accepted_pool);

chain_pool.n_iter = size(Theta_pool,1);
chain_pool.pooled_from = n_chains;
chain_pool.burn_frac_used_for_pooling = burn_frac;
chain_pool.chain_source = chain_source;
chain_pool.source_files = string(chain_files(:));

S = Theta_pool;

chain_pool.posterior_summary = table(chain_pool.param_names(:), ...
    mean(S,1)', ...
    median(S,1)', ...
    prctile(S,2.5,1)', ...
    prctile(S,97.5,1)', ...
    exp(median(S,1))', ...
    'VariableNames', {'parameter','mean_log','median_log','q025_log','q975_log','median_multiplier'});

% -------------------------------------------------------------------------
% Save if requested
% -------------------------------------------------------------------------
if ~isempty(out_file)
    chain = chain_pool;

    out_folder = fileparts(out_file);
    if ~isempty(out_folder) && ~exist(out_folder, 'dir')
        mkdir(out_folder);
    end

    save(out_file, 'chain');

    fprintf('\nSaved pooled chain:\n%s\n', out_file);
end

fprintf('\nPooled %d chains, keeping %.1f%% of each chain.\n', ...
    n_chains, 100*(1-burn_frac));

disp(chain_pool.posterior_summary)

end

%{
qnc_files = {
    'mcmc_outputs/qnc_refractory_50_54_chain1_mcmc51_vp1001_currentSimulator.mat'
    'mcmc_outputs/qnc_refractory_50_54_chain2_mcmc52_vp1001_currentSimulator.mat'
    'mcmc_outputs/qnc_refractory_50_54_chain3_mcmc53_vp1001_currentSimulator.mat'
    'mcmc_outputs/qnc_refractory_50_54_chain4_mcmc54_vp1001_currentSimulator.mat'
};

chain_qnc_pool = pool_chains_from_files( ...
    qnc_files, ...
    0.5, ...
    'mcmc_outputs/pooled_QNC7D_4chains_currentSimulator.mat');

%% OR:

chain_qnc_pool = pool_chains_from_files( ...
    'mcmc_outputs/qnc_refractory_50_54_chain*_mcmc*_vp1001_currentSimulator.mat', ...
    0.5, ...
    'mcmc_outputs/pooled_QNC7D_4chains_currentSimulator.mat');
%}