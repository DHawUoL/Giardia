function chain_pool = pool_chains(chains, burn_frac)

if nargin < 2 || isempty(burn_frac)
    burn_frac = 0.5;
end

chain_pool = chains{1};

Theta_pool = [];
LogPost_pool = [];
Accepted_pool = [];

for c = 1:numel(chains)
    ch = chains{c};

    burn = floor(burn_frac * size(ch.theta,1)) + 1;

    Theta_pool = [Theta_pool; ch.theta(burn:end,:)];
    LogPost_pool = [LogPost_pool; ch.logpost(burn:end)];
    Accepted_pool = [Accepted_pool; ch.accepted(burn:end)];
end

chain_pool.theta = Theta_pool;
chain_pool.logpost = LogPost_pool;
chain_pool.accepted = Accepted_pool;
chain_pool.accept_rate = mean(Accepted_pool);

chain_pool.n_iter = size(Theta_pool,1);
chain_pool.pooled_from = numel(chains);
chain_pool.burn_frac_used_for_pooling = burn_frac;

S = Theta_pool;

chain_pool.posterior_summary = table(chain_pool.param_names(:), ...
    mean(S,1)', ...
    median(S,1)', ...
    prctile(S,2.5,1)', ...
    prctile(S,97.5,1)', ...
    exp(median(S,1))', ...
    'VariableNames', {'parameter','mean_log','median_log','q025_log','q975_log','median_multiplier'});

end

%{
chains = {c1,c2,c3,c4};
chain_pool = pool_chains(chains, 0.5);

disp(chain_pool.posterior_summary)
%}