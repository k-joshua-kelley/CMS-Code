clear; clc; close all;

rng("shuffle");

ts = string(datetime("now", 'Format', 'yyyy-MM-dd_HH-mm-ss'));
save("samples_Nsteps10_Burn100_Thin2_"+ts+".mat", "ts");

data = load("Data/iso_data_stats.mat");
T    = data.T;
f    = data.f;
Nx   = 2^7;
dr   = data.dr;
obs  = data.mu;
kD   = data.kappa;
lnsx = log(2);
lnsy = log(2);

priors = load("iso_priors.mat");
data = load("iso_MAP_results.mat");

hmc_params = load("Tuned_Smp_Params.mat");

smp = hmcSampler( ...
    @(x) log_post_obj(x, kD, obs, lnsx, lnsy, f, Nx, dr, T, priors.priors), ...
    data.x, CheckGradient=false, StepSize=0.025, NumSteps=10, ...
    MassVector=hmc_params.massVector ...
);

[chain,endpoint,accratio] = drawSamples(smp, Burnin=100, NumSamples=500, ThinSize=2, VerbosityLevel=1);
save("samples_Nsteps10_Burn100_Thin2_"+ts+".mat", "chain", "endpoint", "accratio");

function [nlp, grad] = log_post_obj(x, kD, obs, lnsx, lnsy, f, Nx, dr, T, priors_new)
    if nargout < 2
        nlp = iso_nl_post(x(1:7), x(8:14), x(15), x(16), x(17:23), x(24:30), x(31), x(32), x(33), x(34:37), kD, obs, lnsx, lnsy, f, Nx, dr, T, priors_new);
    else
        [nlp, grad] = iso_nl_post(x(1:7), x(8:14), x(15), x(16), x(17:23), x(24:30), x(31), x(32), x(33), x(34:37), kD, obs, lnsx, lnsy, f, Nx, dr, T, priors_new);
        grad = -grad;
    end
    nlp = -nlp;
end