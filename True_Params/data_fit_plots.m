clear; clc; close all;

T = 293:900;

data_fits = load("true_mat_params.mat");

graphite_k_perp = (3*polyval(data_fits.graphite.k_eff, T) - polyval(data_fits.graphite.k_par, T))/2;

figure; hold on;
plot(T, polyval(data_fits.gold_film.k, T))
plot(T, graphite_k_perp)
plot(T, polyval(data_fits.graphite.k_par, T))
plot(T, polyval(data_fits.ZrO2.k, T))

legend("Gold Film", "Graphite (Perp.)", "Graphite (Par.)", "YSZ", Interpreter="latex")
xlabel("Temperature [K]", Interpreter="latex")
ylabel("Thermal Conductivity [W/m/K]", Interpreter="latex")

figure; hold on;
plot(T, polyval(data_fits.gold_film.c, T).*polyval(data_fits.gold_film.rho, T))
plot(T, polyval(data_fits.graphite.c, T).*polyval(data_fits.graphite.rho, T))
plot(T, polyval(data_fits.ZrO2.c, T).*polyval(data_fits.ZrO2.rho, T))

legend("Gold Film", "Graphite", "YSZ", Interpreter="latex")
xlabel("Temperature [K]", Interpreter="latex")
ylabel("Volumetric Heat [J/cc/K]", Interpreter="latex")

figure; hold on;
plot(T, polyval(data_fits.gold_film.k, T)./(polyval(data_fits.gold_film.c, T).*polyval(data_fits.gold_film.rho, T)))
plot(T, graphite_k_perp./(polyval(data_fits.graphite.c, T).*polyval(data_fits.graphite.rho, T)))
plot(T, polyval(data_fits.graphite.k_par, T)./(polyval(data_fits.graphite.c, T).*polyval(data_fits.graphite.rho, T)))
plot(T, polyval(data_fits.ZrO2.k, T)./(polyval(data_fits.ZrO2.c, T).*polyval(data_fits.ZrO2.rho, T)))

legend("Gold Film", "Graphite (Perp.)", "Graphite (Par.)", "YSZ", Interpreter="latex")
xlabel("Temperature [K]", Interpreter="latex")
ylabel("Thermal Diffusivity [mm$^2$/s]", Interpreter="latex")

graphite_k_eff = comsol_fun(data_fits.graphite.k_eff, "[W/m/K]")
graphite_k_par = comsol_fun(data_fits.graphite.k_par, "[W/m/K]")
graphite_c = comsol_fun(data_fits.graphite.c, "[J/g/K]")
graphite_rho = comsol_fun(data_fits.graphite.rho, "[g/cm^3]")

gold_k = comsol_fun(data_fits.gold_film.k, "[W/m/K]")
gold_c = comsol_fun(data_fits.gold_film.c, "[J/g/K]")
gold_rho = comsol_fun(data_fits.gold_film.rho, "[g/cm^3]")

YSZ_k = comsol_fun(data_fits.ZrO2.k, "[W/m/K]")
YSZ_c = comsol_fun(data_fits.ZrO2.c, "[J/g/K]")
YSZ_rho = comsol_fun(data_fits.ZrO2.rho, "[g/cm^3]")

data_fits.gold_film.alpha
data_fits.ZrO2.alpha
data_fits.graphite.alpha

function cmd = comsol_fun(p, units)
    pow = (length(p)-1:-1:0);
    cmd = strjoin(string(compose('%.4g', p)) + units + "*T^" + pow + "[1/K^" + pow + "]", " + ");
end