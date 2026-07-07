% mcmcLoop.m
% Parallel MCMC sensitivity runner.

clear
clc

p = gcp('nocreate');
if isempty(p)
    parpool('local', 4);
end

n_iter     = 3000;
n_vp       = 150;
t_end_days = 35;

cases = { ...
    'pooled_alb_nitro_mtz5_abz5', ...
    'pooled_alb_nitro_mtz7_abz7', ...
    'pooled_alb_nitro_mtz14_abz7', ...
    'pooled_alb_nitro_tdz2g_abz5', ...
    'pooled_alb_nitro_tdz2g_abz7' ...

    
    };

n_cases = numel(cases);
seeds = 42 + (1:n_cases);

outdir = 'mcmc_sensitivity_outputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

results = cell(n_cases,1);

parfor i = 1:n_cases

    case_name = cases{i};
    seed_i = seeds(i);

    outfile = fullfile(outdir, sprintf('%s_chain.mat', case_name));

    if exist(outfile, 'file')
        R = struct();
        R.case_name = case_name;
        R.seed = seed_i;
        R.status = "skipped";
        R.message = "output file already exists";
        R.outfile = string(outfile);
        results{i} = R;
    else
        results{i} = runOneMCMCcase(case_name, seed_i, n_iter, n_vp, t_end_days, outdir);
    end
end

% Convert results to table
case_name = strings(n_cases,1);
seed = zeros(n_cases,1);
status = strings(n_cases,1);
message = strings(n_cases,1);
outfile = strings(n_cases,1);

for i = 1:n_cases
    R = results{i};
    case_name(i) = string(R.case_name);
    seed(i) = R.seed;
    status(i) = string(R.status);
    message(i) = string(R.message);
    outfile(i) = string(R.outfile);
end

run_summary = table(case_name, seed, status, message, outfile);

disp(run_summary)

summary_file = fullfile(outdir, 'run_summary.csv');
writetable(run_summary, summary_file);

fprintf('\nSaved run summary: %s\n', summary_file);