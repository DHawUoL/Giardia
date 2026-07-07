% RUN_ALB_NITRO_SENSITIVITY_MCMC_PARFEVAL
% One-chain structural sensitivity sweep for pooled Alb + nitroimidazole
% regimen mappings.
%
% Outputs are saved into mcmc_sensitivity_outputs.
%
% This assumes giardia_montecarlo has scenario cases for:
%   pooled_alb_nitro_mtz5_abz5
%   pooled_alb_nitro_mtz7_abz7
%   pooled_alb_nitro_tdz2g_abz5
%   pooled_alb_nitro_tdz2g_abz7

clear
clc

%% ------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------

outdir = 'mcmc_sensitivity_outputs';

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

n_iter     = 3000;
n_vp       = 150;
t_end_days = 35;

vp_seed = 1001;

% One chain per scenario for quick structural sensitivity
mcmc_seed_base = 142;

% Four-parameter Alb + nitro calibration:
%   log_nitro_Emax
%   log_nitro_MIC
%   log_ABZ_Emax
%   log_ABZ_MIC
theta_init = [0.0 0.0 0.0 0.0];

jobs = struct([]);

jobs(1).scenario_name = 'pooled_alb_nitro_mtz5_abz5';
jobs(1).regimen       = 'MTZ5D_ABZ5D';
jobs(1).label         = 'MTZ5D + ABZ5D';

jobs(2).scenario_name = 'pooled_alb_nitro_mtz7_abz7';
jobs(2).regimen       = 'MTZ7D_ABZ7D';
jobs(2).label         = 'MTZ7D + ABZ7D';

jobs(3).scenario_name = 'pooled_alb_nitro_tdz2g_abz5';
jobs(3).regimen       = 'TDZ2G_ABZ5D';
jobs(3).label         = 'TDZ2G + ABZ5D';

jobs(4).scenario_name = 'pooled_alb_nitro_tdz2g_abz7';
jobs(4).regimen       = 'TDZ2G_ABZ7D';
jobs(4).label         = 'TDZ2G + ABZ7D';

for j = 1:numel(jobs)
    jobs(j).mcmc_seed = mcmc_seed_base + j - 1;
    jobs(j).vp_seed   = vp_seed;
    jobs(j).theta_init = theta_init;
end

%% ------------------------------------------------------------------------
% Start pool
% -------------------------------------------------------------------------

p = gcp('nocreate');
if isempty(p)
    p = parpool('local');
end

n_jobs = numel(jobs);
futures = parallel.FevalFuture.empty(n_jobs,0);

fprintf('\nSubmitting %d Alb + nitro sensitivity jobs...\n', n_jobs)

for j = 1:n_jobs
    futures(j) = parfeval( ...
        p, ...
        @run_one_sensitivity_chain, ...
        1, ...
        j, ...
        jobs(j), ...
        n_iter, ...
        n_vp, ...
        t_end_days, ...
        outdir);
end

%% ------------------------------------------------------------------------
% Collect results
% -------------------------------------------------------------------------

results = cell(n_jobs,1);

for k = 1:n_jobs
    [completedIdx, result] = fetchNext(futures);
    results{completedIdx} = result;

    fprintf('\nCompleted future %d:\n', completedIdx)
    disp(result)
end

%% ------------------------------------------------------------------------
% Summary table
% -------------------------------------------------------------------------

job_id = zeros(n_jobs,1);
scenario_name = strings(n_jobs,1);
regimen = strings(n_jobs,1);
label = strings(n_jobs,1);
mcmc_seed = zeros(n_jobs,1);
vp_seed_col = zeros(n_jobs,1);
status = strings(n_jobs,1);
outfile = strings(n_jobs,1);
message = strings(n_jobs,1);

for j = 1:n_jobs
    R = results{j};

    job_id(j) = R.job_id;
    scenario_name(j) = R.scenario_name;
    regimen(j) = R.regimen;
    label(j) = R.label;
    mcmc_seed(j) = R.mcmc_seed;
    vp_seed_col(j) = R.vp_seed;
    status(j) = R.status;
    outfile(j) = R.outfile;
    message(j) = R.message;
end

run_summary = table(job_id, scenario_name, regimen, label, ...
    mcmc_seed, vp_seed_col, status, outfile, message);

disp(run_summary)

summary_file = fullfile(outdir, 'alb_nitro_sensitivity_run_summary.csv');
writetable(run_summary, summary_file);

fprintf('\nSaved run summary:\n%s\n', summary_file)

%% ------------------------------------------------------------------------
% Display posterior summaries for successful runs
% -------------------------------------------------------------------------

ok = status == "success";

for j = 1:n_jobs
    if ok(j)
        S = load(outfile(j), 'chain');

        fprintf('\n============================================================\n')
        fprintf('%s\n', label(j))
        fprintf('%s\n', outfile(j))
        fprintf('============================================================\n')

        if isfield(S.chain, 'posterior_summary')
            disp(S.chain.posterior_summary)
        else
            disp(local_posterior_summary(S.chain))
        end

        if isfield(S.chain, 'posterior_median_fit')
            disp(S.chain.posterior_median_fit)
        end
    end
end

%% =========================================================================
% Local worker
% =========================================================================
function R = run_one_sensitivity_chain(job_id, job, n_iter, n_vp, t_end_days, outdir)

    R = struct();
    R.job_id = job_id;
    R.scenario_name = string(job.scenario_name);
    R.regimen = string(job.regimen);
    R.label = string(job.label);
    R.mcmc_seed = job.mcmc_seed;
    R.vp_seed = job.vp_seed;
    R.status = "started";
    R.outfile = "";
    R.message = "";

    try
        chain = giardia_montecarlo( ...
            n_iter, ...
            n_vp, ...
            t_end_days, ...
            job.mcmc_seed, ...
            job.scenario_name, ...
            job.vp_seed, ...
            job.theta_init);

        outname = fullfile(outdir, ...
            sprintf('%s_onechain_mcmc%d_vp%d_currentSimulator.mat', ...
            job.scenario_name, job.mcmc_seed, job.vp_seed));

        scenario_name = job.scenario_name;
        regimen = job.regimen;
        label = job.label;
        mcmc_seed = job.mcmc_seed;
        vp_seed = job.vp_seed;
        theta_init = job.theta_init;

        save(outname, 'chain', 'scenario_name', 'regimen', 'label', ...
            'mcmc_seed', 'vp_seed', 'theta_init');

        R.status = "success";
        R.outfile = string(outname);
        R.message = "completed";

    catch ME

        errname = fullfile(outdir, ...
            sprintf('%s_onechain_mcmc%d_ERROR.mat', ...
            job.scenario_name, job.mcmc_seed));

        scenario_name = job.scenario_name;
        regimen = job.regimen;
        label = job.label;
        mcmc_seed = job.mcmc_seed;
        vp_seed = job.vp_seed;
        theta_init = job.theta_init;

        save(errname, 'ME', 'scenario_name', 'regimen', 'label', ...
            'mcmc_seed', 'vp_seed', 'theta_init');

        R.status = "failed";
        R.outfile = string(errname);
        R.message = string(ME.message);
    end
end

%% =========================================================================
% Local posterior summary fallback
% =========================================================================
function T = local_posterior_summary(chain)

    S = chain.theta;

    T = table(chain.param_names(:), ...
        mean(S,1)', ...
        median(S,1)', ...
        prctile(S,2.5,1)', ...
        prctile(S,97.5,1)', ...
        exp(median(S,1))', ...
        'VariableNames', {'parameter','mean_log','median_log','q025_log','q975_log','median_multiplier'});
end