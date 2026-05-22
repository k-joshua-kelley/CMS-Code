clear; clc;

GRANTA_data = load("GRANTA_material_data.mat");
Folsom_data_matrix_only = load("FolsomCP2012_data_matrix_only.mat");
Folsom_data_surrogate = load("FolsomCP2012_data_Surrogate_TRISO.mat");
optical_data = load("optical_properties.mat");

gold_film.k = polyfit(GRANTA_data.gold.k(:,1), GRANTA_data.gold.k(:,2)/3, 2); % [W/m/K]
gold_film.c = polyfit(GRANTA_data.gold.c(:,1), GRANTA_data.gold.c(:,2)/1e3, 3); % [J/g/K]
gold_e = polyfit(GRANTA_data.gold.expansion(:,1), GRANTA_data.gold.expansion(:,2), 1);

T = 300:900;

for i = 1:length(T)
    gold_rho(i) = 19.32 * exp(-integral(@(Ti) 3*polyval(gold_e, Ti),293,T(i)));
end

gold_film.rho = polyfit(T,gold_rho,1); % [g/cc]


ZrO2.k = polyfit(GRANTA_data.ZrO2.k(:,1), GRANTA_data.ZrO2.k(:,2), 3); % [W/m/K]
ZrO2.c = polyfit(GRANTA_data.ZrO2.c(:,1), GRANTA_data.ZrO2.c(:,2)/1e3, 4); % [J/g/K]
ZrO2_e = polyfit(GRANTA_data.ZrO2.expansion(:,1), GRANTA_data.ZrO2.expansion(:,2), 2);

for i = 1:length(T)
    ZrO2_rho(i) = 6.05 * exp(-integral(@(Ti) 3*polyval(ZrO2_e, Ti),293,T(i)));
end

ZrO2.rho = polyfit(T,ZrO2_rho,1); % [g/cc]


Folsom_data_matrix_only_stacked = [];
for i = 1:length(Folsom_data_matrix_only.matrix_only_data)
    Folsom_data_matrix_only_stacked = [Folsom_data_matrix_only_stacked; Folsom_data_matrix_only.matrix_only_data{i}(:,1) + 273.15, Folsom_data_matrix_only.matrix_only_data{i}(:,2)];
end

graphite.k_par = polyfit(Folsom_data_matrix_only_stacked(:,1), Folsom_data_matrix_only_stacked(:,2), 2);

Folsom_data_surrogate_stacked = [];
for i = 1:length(Folsom_data_surrogate.data_13015)
    Folsom_data_surrogate_stacked = [Folsom_data_surrogate_stacked; Folsom_data_surrogate.data_13015{i}(:,1) + 273.15, Folsom_data_surrogate.data_13015{i}(:,2)];
end

graphite.k_eff = polyfit(Folsom_data_surrogate_stacked(:,1), Folsom_data_surrogate_stacked(:,2), 2);

graphite.c = polyfit(GRANTA_data.graphite.c(:,1), GRANTA_data.graphite.c(:,2)/1e3, 3); % [J/g/K]
graph_e = polyfit(GRANTA_data.graphite.expansion(:,1), GRANTA_data.graphite.expansion(:,2), 1);
for i = 1:length(T)
    graphite_rho(i) = 1.7 * exp(-integral(@(Ti) 3*polyval(graph_e, Ti),293,T(i)));
end

graphite.rho = polyfit(T, graphite_rho, 2);

gold_film.alpha = optical_data.gold.alpha;
ZrO2.alpha = optical_data.ZrO2.alpha;
graphite.alpha = optical_data.graphite.alpha;

save("true_mat_params.mat", "graphite", "gold_film", "ZrO2")

function R2 = r_sq(x,y,d)
p = polyfit(x,y,d);
y_fit = polyval(p, x);       % predicted values

SS_res = sum((y - y_fit).^2);          % residual sum of squares
SS_tot = sum((y - mean(y)).^2);        % total sum of squares

R2 = 1 - SS_res / SS_tot;
end