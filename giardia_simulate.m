function out = giardia_simulate(tvec, regimen, p, opts)
% GIARDIA_SIMULATE  Simulate within-host Giardia dynamics with drug effects.
%
% tvec:    time grid (hours)
% regimen: struct from giardia_regimens
% p:       parameters from giardia_params
% opts:    struct with optional fields:
%          .useMicrobiome       true/false
%          .resistant           true/false
%          .absorbExtinction    true/false
%          .extinction_threshold
%          .extinction_hold_h
%
% Returns struct out with fields:
%   T, C, M, Cpk, clearance_time_h, extinction_time_h, rebounded, regimen_name

if nargin < 4 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'useMicrobiome')
    opts.useMicrobiome = true;
end

if ~isfield(opts,'resistant')
    opts.resistant = false;
end

% -------------------------------------------------------------------------
% Absorbing-extinction switch
% -------------------------------------------------------------------------
% If true, once T has remained below an extinction threshold for long enough,
% the infection is treated as biologically extinct and T=C=0 thereafter.
%
% This is distinct from operational clearance:
%   clearance_threshold: clinical/model detection threshold
%   extinction_threshold: deeper biological extinction threshold
% -------------------------------------------------------------------------
if ~isfield(opts,'absorbExtinction')
    opts.absorbExtinction = false;
end

if ~isfield(opts,'extinction_threshold')
    if isfield(p,'extinction_threshold')
        opts.extinction_threshold = p.extinction_threshold;
    else
        opts.extinction_threshold = 1e-6;
    end
end

if ~isfield(opts,'extinction_hold_h')
    if isfield(p,'extinction_hold_h')
        opts.extinction_hold_h = p.extinction_hold_h;
    else
        opts.extinction_hold_h = 72;
    end
end

Cpk = giardia_pk(tvec, regimen, p);

% Precompute total kill rate over time.
kdrug = zeros(numel(tvec),1);
drugs = fieldnames(Cpk);

for d = 1:numel(drugs)
    nm = drugs{d};

    Emax = p.(nm).Emax;
    if opts.resistant
        Emax = Emax * p.resistance.scale;
    end

    MIC = p.(nm).MIC;
    h   = p.(nm).h;
    Ce  = Cpk.(nm).eff;

    kdrug = kdrug + Emax .* (Ce.^h) ./ (MIC.^h + Ce.^h);
end

% ODE integration: forward Euler.
dt = tvec(2) - tvec(1);

T = zeros(size(tvec));
C = zeros(size(tvec));
M = zeros(size(tvec));

T(1) = p.T0;
C(1) = p.C0;
M(1) = p.M0;

extinct = false;
extinction_time_h = NaN;
below_extinction_start_h = NaN;

for i = 1:numel(tvec)-1

    % Once absorbed, remain extinct.
    if extinct
        T(i+1) = 0;
        C(i+1) = 0;

        % Microbiome can still evolve after parasite extinction.
        dM = p.rM*M(i)*(1 - M(i)) - p.dys*M(i);
        M(i+1) = min(1, max(0, M(i) + dt*dM));

        continue
    end

    km = p.k_micro;
    if ~opts.useMicrobiome
        km = 0;
    end

    % Simple ramp for immune clearance after 7 days of infection history.
    % tvec may begin before treatment, so use time since simulation start
    % rather than absolute model time.
    infection_age_h = tvec(i) - tvec(1) + 24*30;
    kI = p.k_I * (1 + 0.5*(infection_age_h > 24*7));

    dT = p.rT*T(i)*(1 - T(i)/p.K) ...
       - p.k_encyst*T(i) ...
       - kI*T(i) ...
       - km*M(i)*T(i) ...
       - kdrug(i)*T(i);

    dC = p.k_encyst*T(i) - p.k_c_decay*C(i);

    dM = p.rM*M(i)*(1 - M(i)) - p.dys*M(i);

    T(i+1) = max(0, T(i) + dt*dT);
    C(i+1) = max(0, C(i) + dt*dC);
    M(i+1) = min(1, max(0, M(i) + dt*dM));

    % ---------------------------------------------------------------------
    % Absorbing-extinction rule.
    % ---------------------------------------------------------------------
    if opts.absorbExtinction

        if T(i+1) < opts.extinction_threshold

            if isnan(below_extinction_start_h)
                below_extinction_start_h = tvec(i+1);
            end

            if (tvec(i+1) - below_extinction_start_h) >= opts.extinction_hold_h
                extinct = true;
                extinction_time_h = below_extinction_start_h;

                T(i+1) = 0;
                C(i+1) = 0;
            end

        else
            below_extinction_start_h = NaN;
        end
    end
end

% Clearance time: first POST-TREATMENT time T<thr and stays below
% for the required hold period. Treatment begins at t = 0.
thr = p.clearance_threshold;
hold = p.clearance_hold_h;
clear_time = NaN;

idx_treat = find(tvec >= 0, 1, 'first');

for i = idx_treat:numel(tvec)
    if T(i) < thr
        j = i + round(hold/dt);

        if j <= numel(tvec) && all(T(i:j) < thr)
            clear_time = tvec(i);
            break
        end
    end
end

% Rebound flag: after first operational clearance, does T ever rise above threshold again?
rebounded = false;

if ~isnan(clear_time)
    idx_clear = find(tvec >= clear_time, 1, 'first');

    if ~isempty(idx_clear)
        rebounded = any(T(idx_clear:end) > thr);
    end
end

out.T = T;
out.C = C;
out.M = M;
out.Cpk = Cpk;
out.clearance_time_h = clear_time;
out.extinction_time_h = extinction_time_h;
out.rebounded = rebounded;
out.regimen_name = regimen.name;
out.tvec = tvec(:);
out.kdrug = kdrug;
out.opts = opts;

end