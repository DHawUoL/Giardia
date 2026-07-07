% RUN_QNC_PRIOR_SENSITIVITY_MCMC_PARFEVAL
% One-chain QNC prior-width sensitivity sweep.
%
% Requires giardia_montecarlo to accept:
%   prior_sd_override as optional 8th argument.

clear
clc

outdir = 'mcmc_sensitivity_outputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

scenario_name = 'qnc_refractory_50_54';

n_iter     = 3000;
n_vp       = 150;
t_end_days = 35;

vp_seed = 1001;

theta_init = [0.0 0.0];

jobs = struct([]);

jobs(1).prior_sd = [0.60 0.60];
jobs(1).label = 'QNC prior sd 0.60';

jobs(2).prior_sd = [0.90 0.90];
jobs(2).label = 'QNC prior sd 0.90 reference';

jobs(3).prior_sd = [1.25 1.25];
jobs(3).label = 'QNC prior sd 1.25';

jobs(4).prior_sd = [1.50 1.50];
jobs(4).label = 'QNC prior sd 1.50';

mcmc_seed_base = 251;

for j = 1:numel(jobs)
    jobs(j).scenario_name = scenario_name;
    jobs(j).regimen = 'QNC7D';
    jobs(j).mcmc_seed = mcmc_seed_base + j - 1;
    jobs(j).vp_seed = vp_seed;
    jobs(j).theta_init = theta_init;
end

p = gcp('nocreate');
if isempty(p)
    p = parpool('local');
end

n_jobs = numel(jobs);
futures = parallel.FevalFuture.empty(n_jobs,0);

fprintf('\nSubmitting %d QNC prior-sensitivity jobs...\n', n_jobs)

for j = 1:n_jobs
    futures(j) = parfeval( ...
        p, ...
        @run_one_qnc_prior_chain, ...
        1, ...
        j, ...
        jobs(j), ...
        n_iter, ...
        n_vp, ...
        t_end_days, ...
        outdir);
end

results = cell(n_jobs,1);

for k = 1:n_jobs
    [completedIdx, result] = fetchNext(futures);
    results{completedIdx} = result;

    fprintf('\nCompleted future %d:\n', completedIdx)
    disp(result)
end

job_id = zeros(n_jobs,1);
label = strings(n_jobs,1);
prior_sd = strings(n_jobs,1);
mcmc_seed = zeros(n_jobs,1);
vp_seed_col = zeros(n_jobs,1);
status = strings(n_jobs,1);
outfile = strings(n_jobs,1);
message = strings(n_jobs,1);

for j = 1:n_jobs
    R = results{j};

    job_id(j) = R.job_id;
    label(j) = R.label;
    prior_sd(j) = mat2str(R.prior_sd);
    mcmc_seed(j) = R.mcmc_seed;
    vp_seed_col(j) = R.vp_seed;
    status(j) = R.status;
    outfile(j) = R.outfile;
    message(j) = R.message;
end

run_summary = table(job_id, label, prior_sd, mcmc_seed, vp_seed_col, ...
    status, outfile, message);

disp(run_summary)

summary_file = fullfile(outdir, 'qnc_prior_sensitivity_run_summary.csv');
writetable(run_summary, summary_file);

fprintf('\nSaved run summary:\n%s\n', summary_file)

% Display posterior summaries
for j = 1:n_jobs
    if status(j) == "success"
        S = load(outfile(j), 'chain');

        fprintf('\n============================================================\n')
        fprintf('%s\n', label(j))
        fprintf('%s\n', outfile(j))
        fprintf('============================================================\n')

        disp(S.chain.posterior_summary)

        if isfield(S.chain, 'posterior_median_fit')
            local_disp_fit(S.chain.posterior_median_fit)
        end
    end
end

%% =========================================================================
function R = run_one_qnc_prior_chain(job_id, job, n_iter, n_vp, t_end_days, outdir)

    R = struct();
    R.job_id = job_id;
    R.label = string(job.label);
    R.scenario_name = string(job.scenario_name);
    R.regimen = string(job.regimen);
    R.prior_sd = job.prior_sd;
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
            job.theta_init, ...
            job.prior_sd);

        prior_tag = sprintf('sd%.2f_%.2f', job.prior_sd(1), job.prior_sd(2));
        prior_tag = strrep(prior_tag, '.', 'p');

        outname = fullfile(outdir, ...
            sprintf('qnc_refractory_50_54_prior_%s_onechain_mcmc%d_vp%d_currentSimulator.mat', ...
            prior_tag, job.mcmc_seed, job.vp_seed));

        scenario_name = job.scenario_name;
        regimen = job.regimen;
        label = job.label;
        prior_sd_override = job.prior_sd;
        mcmc_seed = job.mcmc_seed;
        vp_seed = job.vp_seed;
        theta_init = job.theta_init;

        save(outname, 'chain', 'scenario_name', 'regimen', 'label', ...
            'prior_sd_override', 'mcmc_seed', 'vp_seed', 'theta_init');

        R.status = "success";
        R.outfile = string(outname);
        R.message = "completed";

    catch ME

        prior_tag = sprintf('sd%.2f_%.2f', job.prior_sd(1), job.prior_sd(2));
        prior_tag = strrep(prior_tag, '.', 'p');

        errname = fullfile(outdir, ...
            sprintf('qnc_refractory_50_54_prior_%s_onechain_mcmc%d_ERROR.mat', ...
            prior_tag, job.mcmc_seed));

        scenario_name = job.scenario_name;
        regimen = job.regimen;
        label = job.label;
        prior_sd_override = job.prior_sd;
        mcmc_seed = job.mcmc_seed;
        vp_seed = job.vp_seed;
        theta_init = job.theta_init;

        save(errname, 'ME', 'scenario_name', 'regimen', 'label', ...
            'prior_sd_override', 'mcmc_seed', 'vp_seed', 'theta_init');

        R.status = "failed";
        R.outfile = string(errname);
        R.message = string(ME.message);
    end
end

%% =========================================================================
function local_disp_fit(fit_detail)

    if isfield(fit_detail, 'target_rows')
        tr = fit_detail.target_rows;
        if isstruct(tr)
            disp(struct2table(tr))
        else
            disp(tr)
        end
    else
        disp(fit_detail)
    end
end