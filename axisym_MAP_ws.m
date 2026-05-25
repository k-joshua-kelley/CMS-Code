clear; clc; close all;

rng("shuffle");
seed = rng();

data = load("Data/axisym_data_stats.mat");
T = data.T;

material_fits = load("True_Params/true_mat_params.mat");
lnkf = log(polyval(material_fits.gold_film.k, T));
lnCf = log(polyval(material_fits.gold_film.c, T).*polyval(material_fits.gold_film.rho, T));
lnaf = log(material_fits.gold_film.alpha);
lnhf = log(0.15);
lnks_perp = log((3*polyval(material_fits.graphite.k_eff, T) - polyval(material_fits.graphite.k_par, T))/2);
lnks_par = log(polyval(material_fits.graphite.k_par, T));
lnCs = log(polyval(material_fits.graphite.c, T) .* polyval(material_fits.graphite.rho, T));
lnas = log(material_fits.graphite.alpha);
lnRth = -inf;

Os = load("True_Params/gold_graphite_O.mat");
Os = Os.v;

lnsx = log(2);
lnsy = log(2);
f = data.f;
Nx = 160;
X_probe = data.X_probe;

tic
Df = exp(lnkf-lnCf); % mm^2/s
Ds = exp(0.5*(lnks_perp+lnks_par)-lnCs); % mm^2/s
Lthf = sqrt(Df ./ pi ./ f); % um
Lths = sqrt(Ds ./ pi ./ f); % um
x_max = max(max(sqrt(X_probe(:,1).^2+X_probe(:,2).^2))*10, 1*max(Lthf, Lths));

load("axisym_priors.mat")

nl_post_handle = @(x) axisym_nl_post(x(1:7),x(8:14),lnaf,x(15),x(16:22),x(23:29),x(30:36),lnas,lnRth,x(37:46),x(47:56),x(57:66),lnsx,lnsy,f,Nx,X_probe,x_max, data.mu, x(67), data.kappa, x(68:72), T, priors);

nl_post_warm_start_handle = @(x,i) axisym_nl_post_warm_start(x(1),x(2),lnaf,x(3),x(4),x(5),x(6),lnas,lnRth,x(7),x(8),x(9),lnsx,lnsy,f,Nx,X_probe,x_max, data.mu(i,:,:,:), 1e6, data.kappa(i,:,:,:), priors);

options = optimoptions("fmincon", "FiniteDifferenceType","central", Display="iter-detailed",SpecifyConstraintGradient=true, SpecifyObjectiveGradient=true);

x_ws = cell(length(T),1);
fval_ws = cell(length(T),1);
exitflag_ws = cell(length(T),1);
output_ws = cell(length(T),1);
lambda_ws = cell(length(T),1);
grad_ws = cell(length(T),1);

Os_theta = normrnd(0, deg2rad(20));
Os_phi = pi*rand;
Os01 = sin(Os_theta) .* cos(Os_phi);
Os02 = cos(Os_theta);
Os03 = sin(Os_theta) .* sin(Os_phi);
for i = 1:length(T)
    x0_warm_start = [
        normrnd(priors.lnkf.mu, priors.lnkf.sigma), ...
        normrnd(priors.lnCf.mu, priors.lnCf.sigma), ...
        normrnd(priors.lnhf.mu, priors.lnhf.sigma), ...
        normrnd(priors.lnks_perp.mu, priors.lnks_perp.sigma), ...
        normrnd(priors.lnks_par.mu, priors.lnks_par.sigma), ...
        normrnd(priors.lnCs.mu, priors.lnCs.sigma), ...
        Os01, Os02, Os03
    ];
    [x_ws{i},fval_ws{i},exitflag_ws{i},output_ws{i},lambda_ws{i},grad_ws{i}] = fmincon(@(x) nl_post_warm_start_handle(x,i), x0_warm_start, [], [], [], [], [], [], @nonlcon_warm_start, options);
    x_ws{i} = x_ws{i}(:);
end
x_ws = horzcat(x_ws{:}).';
save("axisym_MAP_results_ws"+seed.Seed+".mat", "x_ws", "fval_ws", "exitflag_ws", "output_ws", "lambda_ws", "grad_ws", "seed");

function [ineq, eq, Gineq, Geq] = nonlcon_warm_start(x)
arguments
    x (:,1)
end

    ineq = [];
    Gineq = [];

    eq = x(7:9).' * x(7:9) - 1;

    Geq = zeros(length(x),1);
    Geq(7:9) = 2*x(7:9);

end