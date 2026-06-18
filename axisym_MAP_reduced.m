clear; clc; close all;

load("axisym_MAP_results1660744828.mat")
load("axisym_MAP_results_reduced(6)1660744828.mat")

rng(seed2.Seed);

data = load("Data/axisym_data_stats.mat");
T = data.T;

material_fits = load("True_Params/true_mat_params.mat");
lnkf = log(polyval(material_fits.gold_film.k, T));
lnCf = log(polyval(material_fits.gold_film.c, T).*polyval(material_fits.gold_film.rho, T));
lnaf = log(material_fits.gold_film.alpha);
lnks_perp = log((3*polyval(material_fits.graphite.k_eff, T) - polyval(material_fits.graphite.k_par, T))/2);
lnks_par = log(polyval(material_fits.graphite.k_par, T));
lnCs = log(polyval(material_fits.graphite.c, T) .* polyval(material_fits.graphite.rho, T));
lnas = log(material_fits.graphite.alpha);
lnRth = -inf;

lnsx = log(2);
lnsy = log(2);
f = data.f;
Nx = 160;
X_probe = data.X_probe;

Df = exp(lnkf-lnCf); % mm^2/s
Ds = exp(0.5*(lnks_perp+lnks_par)-lnCs); % mm^2/s
Lthf = sqrt(Df ./ pi ./ f); % um
Lths = sqrt(Ds ./ pi ./ f); % um
x_max = max(max(sqrt(X_probe(:,1).^2+X_probe(:,2).^2))*10, 1*max(Lthf, Lths));

load("axisym_priors.mat")

Oi = randperm(10,5);

nl_post_handle = @(x) axisym_nl_post_helper(x(1:7),x(8:14),lnaf,x(15),x(16:22),x(23:29),x(30:36),lnas,lnRth,x(37:41),x(42:46),x(47:51),lnsx,lnsy,f,Nx,X_probe,x_max, data.mu(:,Oi,:,:), x(52), data.kappa(:,Oi,:,:), x(53:57), T, priors);

options = optimoptions("fmincon", ...
    FiniteDifferenceType="central", ...
    Display="iter-detailed", ...
    SpecifyConstraintGradient=true, ...
    SpecifyObjectiveGradient=true, ...
    MaxFunctionEvaluations=9000, ...
    Algorithm="interior-point", ...
    EnableFeasibilityMode=true, ...
    SubproblemAlgorithm="cg"...
);

x0 = zeros(57,1);
x0([1:51,53:end]) = x;
x0(52) = log(1e4);

[psi, grad] = nl_post_handle(x0);
%%

[~, err] = checkGradients(@check_nonlcon, x0, options, "Display","on");
max(abs(err.Objective(:)))
[~, err] = checkGradients(nl_post_handle, x0, options, "Display","on");
max(abs(err.Objective(:)))
[x,fval,exitflag,output,lambda,grad] = fmincon(nl_post_handle, x0, [], [], [], [], [], [], @nonlcon, options);

H = central_diff_hessian(nl_post_handle, x, 1e-4);

save("axisym_MAP_results_reduced(7)"+seed2.Seed+".mat", "x", "fval", "exitflag", "output", "lambda", "grad", "H");

function [eq, Jeq] = check_nonlcon(x)
    [~, eq, ~, Geq] = nonlcon(x);
    Jeq = Geq.';
end
function [ineq, eq, Gineq, Geq] = nonlcon(x)

arguments
    x (:,1)
end

ineq = [];
Gineq = [];

% Construct matrix
O = [x(37:41), x(42:46), x(47:51)];

% One equality constraint per row:
% sum(row.^2) - 1 = 0
eq = sum(O.^2, 2) - 1;

if nargout > 3

    % Number of variables and constraints
    nx = length(x);
    neq = length(eq);

    Geq = zeros(nx, neq);

    % Gradients
    %
    % Row i constraint:
    % O(i,1)^2 + O(i,2)^2 + O(i,3)^2 - 1
    %
    % derivative wrt each involved variable = 2*x

    for i = 1:neq
        Geq(36+i, i) = 2*x(36+i); % first column
        Geq(41+i, i) = 2*x(41+i); % second column
        Geq(46+i, i) = 2*x(46+i); % third column
    end

end
end

function [psi, grad] = axisym_nl_post_helper(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, tau, kD, theta, T, priors)
    kT = exp(tau);
    try
        if nargout < 2
            psi = axisym_nl_post(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD, theta, T, priors);
            psi = psi + nln(tau, priors.tau.mu, priors.tau.sigma);
        else
            [psi, grad] = axisym_nl_post(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD, theta, T, priors);
            [psi_tau, grad_tau] = nln(tau, priors.tau.mu, priors.tau.sigma, [true,false,false]);

            psi = psi + psi_tau;

            grad(52) = grad(52) + grad_tau{1};
        end
    catch err
        ts = string(datetime("now", 'Format', 'yyyy-MM-dd_HH-mm-ss'));
        save("error_" + ts + ".mat", "lnkf", "lnCf", "lnaf", "lnhf", "lnks_perp", "lnks_par", "lnCs", "lnas", "lnRth", "Os1", "Os2", "Os3", "lnsx", "lnsy", "f", "Nx", "X_probe", "x_max", "obs", "kT", "kD", "theta", "T", "priors", "err")
        psi = NaN;
        grad = NaN(56,1);
    end
end