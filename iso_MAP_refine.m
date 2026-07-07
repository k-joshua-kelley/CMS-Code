clear; clc; close all;

% Experimental Parameters (known)
data = load("Data/iso_data_stats.mat");
T    = data.T;
f    = data.f;
Nx   = 160;
dr   = data.dr;
obs  = data.mu;
kD   = 100*data.kappa;
lnsx = log(2);
lnsy = log(2);

true_params = load("True_Params\optical_properties.mat");
lnaf = log(true_params.gold.alpha);
lnas = log(true_params.ZrO2.alpha);
lnRth = -inf;

% Load Priors
load("iso_priors.mat");

% format objective functions for Temp-coupled and Temp-independent
% objectives
nl_post_obj = @(x) iso_nl_post(x(1:7), x(8:14), lnaf, x(15), x(16:22), x(23:29), lnas, lnRth, x(30), x(31:34), kD, obs, lnsx, lnsy, f, Nx*2, dr, T, priors);

options = optimoptions("fminunc", Display="iter-detailed", SpecifyObjectiveGradient=true, FiniteDifferenceType="central", StepTolerance=1e-10);

%% Temperature-Coupled
MAP = load("iso_MAP_results.mat");
x0 = MAP.x;

[~, err] = checkGradients(nl_post_obj, x0, options, "Display","on");
[x,fval,exitflag,output,grad,hessian] = fminunc(nl_post_obj, x0, options);

save("iso_MAP_results_refined.mat", "x", "fval", "exitflag", "output", "grad", "hessian")

function [psi, grad] = nl_post_ind(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, tau, kD, obs, lnsx, lnsy, f, Nx, dr, priors)
    kT = exp(tau);
    % random variables: lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth,
    % kT, theta

    % disp(compose("x = [%s, %s, %s, %s, %s, %s, %s, %s, %s, %s]", ...
    %     join(string(lnkf), ", "), ...
    %     join(string(lnCf), ", "), ...
    %     join(string(lnaf), ", "), ...
    %     join(string(lnhf), ", "), ...
    %     join(string(lnks), ", "), ...
    %     join(string(lnCs), ", "), ...
    %     join(string(lnas), ", "), ...
    %     join(string(lnRth), ", "), ...
    %     join(string(kT), ", "), ...
    %     join(string(theta), ", ")))
    if nargout < 2
        psi_nll      = iso_nll_M(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, kT, kD, obs, lnsx, lnsy, f, Nx, dr);
        phi_GP_lnkf  = nln(lnkf, priors.lnkf.mu, priors.lnkf.sigma, [true, false, false]);
        phi_GP_lnCf  = nln(lnCf, priors.lnCf.mu, priors.lnCf.sigma, [true, false, false]);
        phi_GP_lnhf  = nln(lnhf, priors.lnhf.mu, priors.lnhf.sigma, [true, false, false]);
        phi_GP_lnks  = nln(lnks, priors.lnks.mu, priors.lnks.sigma, [true, false, false]);
        phi_GP_lnCs  = nln(lnCs, priors.lnCs.mu, priors.lnCs.sigma, [true, false, false]);
        psi_tau = nln(tau, priors.tau.mu, priors.tau.sigma, [true,false,false]);
    else
        [psi_nll, grad_nll] = iso_nll_M(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, kT, kD, obs, lnsx, lnsy, f, Nx, dr);
        [phi_GP_lnkf, grad_GP_lnkf]   = nln(lnkf, priors.lnkf.mu, priors.lnkf.sigma, [true, false, false]);
        [phi_GP_lnCf, grad_GP_lnCf]   = nln(lnCf, priors.lnCf.mu, priors.lnCf.sigma, [true, false, false]);
        [phi_GP_lnhf, grad_GP_lnhf]   = nln(lnhf, priors.lnhf.mu, priors.lnhf.sigma, [true, false, false]);
        [phi_GP_lnks, grad_GP_lnks]   = nln(lnks, priors.lnks.mu, priors.lnks.sigma, [true, false, false]);
        [phi_GP_lnCs, grad_GP_lnCs]   = nln(lnCs, priors.lnCs.mu, priors.lnCs.sigma, [true, false, false]);
        [psi_tau, grad_tau] = nln(tau, priors.tau.mu, priors.tau.sigma, [true,false,false]);

        grad_GP = [grad_GP_lnkf{1}; grad_GP_lnCf{1}; 0; grad_GP_lnhf{1}; grad_GP_lnks{1}; grad_GP_lnCs{1}; 0; 0; 0];

        grad = grad_nll + grad_GP;
        grad(9) = kT*grad(9) + grad_tau{1};
        grad = grad([1:2,4:6,9]);
        
        % if any(isnan(grad))
        %     pause(1);
        % end
    end
    psi = psi_nll + phi_GP_lnkf + phi_GP_lnCf + phi_GP_lnhf + phi_GP_lnks + phi_GP_lnCs + psi_tau;
    % if isnan(psi)
    %     pause(1);
    % end
end