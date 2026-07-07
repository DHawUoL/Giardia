%function f = run_mcmc_parallel
% Parallel QNC MCMC chains using parfeval rather than parfor.
% This avoids MATLAB parfor transparency problems.

clear
clc

%{
scenario_name = 'qnc_refractory_50_54';
pooled_name = 'pooled_QNC7D_4chains_currentSimulator.mat';

mcmc_seeds = [51 52 53 54];

theta_inits = [
     0.0   0.0
     1.0   0.0
    -1.0   0.0
     0.5  -0.5
];
%}

scenario_name = 'pooled_alb_nitro_mtz14_abz7';
pooled_name = 'pooled_MTZ14D_ABZ7D_4chains_currentSimulator.mat';

mcmc_seeds = [42 43 44 45];

theta_inits = [
     0.0   0.0   0.0   0.0
     1.0   0.0   1.0   0.0
    -1.0   0.0  -1.0   0.0
     0.5   0.5  -0.5  -0.5
];

n_iter     = 3000;
n_vp       = 150;
t_end_days = 35;

vp_seed = 1001;

outdir = 'mcmc_outputs';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Start pool if needed
p = gcp('nocreate');
if isempty(p)
    p = parpool('local');
end

n_chains = numel(mcmc_seeds);
futures = parallel.FevalFuture.empty(n_chains,0);

fprintf('\nSubmitting %d chains for scenario: %s\n', n_chains, scenario_name)

for c = 1:n_chains
    futures(c) = parfeval( ...
        p, ...
        @run_one_qnc_chain, ...
        1, ...
        c, ...
        scenario_name, ...
        n_iter, ...
        n_vp, ...
        t_end_days, ...
        mcmc_seeds(c), ...
        vp_seed, ...
        theta_inits(c,:), ...
        outdir);
end

results = cell(n_chains,1);

for k = 1:n_chains
    [completedIdx, result] = fetchNext(futures);
    results{completedIdx} = result;

    fprintf('\nCompleted future %d:\n', completedIdx)
    disp(result)
end

% Convert results to table
chain_id = zeros(n_chains,1);
mcmc_seed_col = zeros(n_chains,1);
vp_seed_col = zeros(n_chains,1);
status = strings(n_chains,1);
outfile = strings(n_chains,1);
message = strings(n_chains,1);

for c = 1:n_chains
    R = results{c};

    chain_id(c) = R.chain_id;
    mcmc_seed_col(c) = R.mcmc_seed;
    vp_seed_col(c) = R.vp_seed;
    status(c) = R.status;
    outfile(c) = R.outfile;
    message(c) = R.message;
end

run_summary = table(chain_id, mcmc_seed_col, vp_seed_col, ...
    status, outfile, message);

disp(run_summary)

writetable(run_summary, fullfile(outdir, ...
    sprintf('%s_parfeval_run_summary.csv', scenario_name)));

% Pool successful chains
ok = status == "success";

if sum(ok) >= 2

    chains_qnc = cell(sum(ok),1);
    labels_qnc = cell(sum(ok),1);

    jj = 0;
    for c = 1:n_chains
        if ok(c)
            jj = jj + 1;

            S = load(outfile(c), 'chain');
            chains_qnc{jj} = S.chain;
            labels_qnc{jj} = sprintf('chain %d, seed %d', ...
                chain_id(c), mcmc_seed_col(c));
        end
    end

    chain_qnc_pool = pool_chains(chains_qnc, 0.5);

    chain = chain_qnc_pool;
    pooled_file = fullfile(outdir, pooled_name);

    save(pooled_file, 'chain');

    fprintf('\nSaved pooled QNC chain:\n%s\n', pooled_file)
    disp(chain.posterior_summary)

else
    warning('Fewer than two successful chains. Not pooling.')
end


% -------------------------------------------------------------------------
% Local worker function
% -------------------------------------------------------------------------
function R = run_one_qnc_chain(chain_id, scenario_name, n_iter, n_vp, ...
    t_end_days, mcmc_seed, vp_seed, theta_init, outdir)

    R = struct();
    R.chain_id = chain_id;
    R.scenario_name = scenario_name;
    R.mcmc_seed = mcmc_seed;
    R.vp_seed = vp_seed;
    R.theta_init = theta_init;
    R.status = "started";
    R.outfile = "";
    R.message = "";

    try
        chain = giardia_montecarlo( ...
            n_iter, ...
            n_vp, ...
            t_end_days, ...
            mcmc_seed, ...
            scenario_name, ...
            vp_seed, ...
            theta_init);

        % Use our own clear filename as well as whatever giardia_montecarlo saves.
        outname = fullfile(outdir, ...
            sprintf('%s_chain%d_mcmc%d_vp%d_currentSimulator.mat', ...
            scenario_name, chain_id, mcmc_seed, vp_seed));

        save(outname, 'chain', 'scenario_name', 'chain_id', ...
            'mcmc_seed', 'vp_seed', 'theta_init');

        R.status = "success";
        R.outfile = string(outname);
        R.message = "completed";

    catch ME
        errname = fullfile(outdir, ...
            sprintf('qnc_refractory_50_54_chain%d_mcmc%d_ERROR.mat', ...
            chain_id, mcmc_seed));

        save(errname, 'ME', 'scenario_name', 'chain_id', ...
            'mcmc_seed', 'vp_seed', 'theta_init');

        R.status = "failed";
        R.outfile = string(errname);
        R.message = string(ME.message);
    end
end