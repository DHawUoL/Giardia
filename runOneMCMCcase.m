function result = runOneMCMCcase(case_name, seed_i, n_iter, n_vp, t_end_days, outdir)
% RUNONEMCMCCASE  Run one MCMC sensitivity case and save result.
%
% This is deliberately a separate function so that parfor does not complain
% about transparency violations from save / try / catch logic in the script.

result = struct();
result.case_name = case_name;
result.seed = seed_i;
result.status = "started";
result.message = "";
result.outfile = "";

outfile = fullfile(outdir, sprintf('%s_chain.mat', case_name));
errfile = fullfile(outdir, sprintf('%s_ERROR.mat', case_name));

try
    fprintf('\nStarting: %s\n', case_name);

    chain = giardia_montecarlo( ...
        n_iter, ...
        n_vp, ...
        t_end_days, ...
        seed_i, ...
        case_name);

    S = struct();
    S.chain = chain;
    S.case_name = case_name;
    S.seed = seed_i;

    save(outfile, '-struct', 'S');

    result.status = "success";
    result.message = "completed";
    result.outfile = string(outfile);

catch ME
    E = struct();
    E.message = ME.message;
    E.identifier = ME.identifier;
    E.stack = ME.stack;
    E.case_name = case_name;
    E.seed = seed_i;

    save(errfile, '-struct', 'E');

    result.status = "failed";
    result.message = string(ME.message);
    result.outfile = string(errfile);
end

end