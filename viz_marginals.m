%% Main GP Results
clear; clc; close all;
iso_results = load("iso_MAP_results.mat");
iso_priors = load("iso_priors.mat");
true_params = load("True_Params\true_mat_params.mat");
axisym_results = load("axisym_laplace_params.mat");
axisym_priors = load("axisym_priors.mat");
colors = colororder;
color_true = colors(2,:);
color_post = colors(1,:);
color_samp = [0.5,0.5,0.5];
lw = 1.5; cs = 5; ms = 2.5;

T = (300:100:900).';
Ts = (290:910).';
N_samp = 10;

fprintf("ISO TABLE:\n\n")
for i = [15,30:34]
    indx = false(size(iso_results.hessian,1),1);
    indx(i) = true;
    Hp = schur_comp(iso_results.hessian, indx);
    mode = exp(iso_results.x(i)-1/Hp);
    S = get_logN_HDI_param(sigma=sqrt(1/Hp), alpha=0.95);
    fprintf("$\\mathrm{logN}(\\num{%.4g},\\num{%.4g})$ & \\num{%.4g} & [\\num{%.4g},\\num{%.4g}] \\\\\n", iso_results.x(i), 1/Hp, mode, mode/S, mode*S)
end
disp("\hline")
fprintf("\n\n")


fprintf("AXISYM TABLE:\n\n")
i = 15;
indx = false(size(axisym_results.HessM,1),1);
indx(i) = true;
Hp = schur_comp(axisym_results.HessM, indx);
mode = exp(axisym_results.x(i)-1/Hp);
S = get_logN_HDI_param(sigma=sqrt(1/Hp), alpha=0.95);
fprintf("$\\mathrm{logN}(\\num{%.4g},\\num{%.4g})$ & \\num{%.4g} & [\\num{%.4g},\\num{%.4g}] \\\\\n", axisym_results.x(i), 1/Hp, mode, mode/S, mode*S)
for i = (52:57)-5
    indx = false(size(axisym_results.HessM,1),1);
    indx(i) = true;
    Hp = schur_comp(axisym_results.HessM, indx);
    mode = exp(axisym_results.x(i+5)-1/Hp);
    S = get_logN_HDI_param(sigma=sqrt(1/Hp), alpha=0.95);
    fprintf("$\\mathrm{logN}(\\num{%.4g},\\num{%.4g})$ & \\num{%.4g} & [\\num{%.4g},\\num{%.4g}] \\\\\n", axisym_results.x(i+5), 1/Hp, mode, mode/S, mode*S)
end
disp("\hline")
fprintf("\n\n")

% Iso Samples
indx = false(size(iso_results.hessian,1),1);
indx([1:14,16:29,31:34]) = true;
H_marg = schur_comp(iso_results.hessian, indx);
Sp = inv(H_marg);
Sp = 0.5*(Sp+Sp.');
mu = iso_results.x(indx);
mu0 = @(T) [iso_priors.priors.lnkf.mu, iso_priors.priors.lnCf.mu, iso_priors.priors.lnks.mu, iso_priors.priors.lnCs.mu] .* ones(length(T),1);
sigma0 = {...
    @(T) iso_priors.priors.lnkf.sigma .* ones(length(T),1);...
    @(T) iso_priors.priors.lnCf.sigma .* ones(length(T),1);...
    @(T) iso_priors.priors.lnks.sigma .* ones(length(T),1);...
    @(T) iso_priors.priors.lnCs.sigma .* ones(length(T),1);...
};
iso_samples = sample_posterior_predictive(T, Ts, mu0, sigma0, mu, Sp, N_samp);
iso_samples = reshape(iso_samples, N_samp, [], 4);


% Axisym Samples
indx = false(size(axisym_results.HessM,1),1);
indx([1:14, 16:36, (53:57)-5]) = true;
H_marg = schur_comp(axisym_results.HessM, indx);
Sp = inv(H_marg);
Sp = 0.5*(Sp+Sp.');
mu = axisym_results.x([1:14, 16:36, 53:57]);
mu0 = @(T) [axisym_priors.priors.lnkf.mu, axisym_priors.priors.lnCf.mu, axisym_priors.priors.lnks_perp.mu, axisym_priors.priors.lnks_par.mu, axisym_priors.priors.lnCs.mu] .* ones(length(T),1);
sigma0 = {...
    @(T) axisym_priors.priors.lnkf.sigma .* ones(length(T),1);...
    @(T) axisym_priors.priors.lnCf.sigma .* ones(length(T),1);...
    @(T) axisym_priors.priors.lnks_perp.sigma .* ones(length(T),1);...
    @(T) axisym_priors.priors.lnks_par.sigma .* ones(length(T),1);...
    @(T) axisym_priors.priors.lnCs.sigma .* ones(length(T),1);...
};
axisym_samples = sample_posterior_predictive(T, Ts, mu0, sigma0, mu, Sp, N_samp);
axisym_samples = reshape(axisym_samples, N_samp, [], 5);


figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');



indx = false(size(iso_results.hessian,1),1);
indx(1:14) = true;
Hp = schur_comp(iso_results.hessian, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');
Sp = Sp(1:7,1:7) + Sp(8:14,8:14) - 2*Sp(1:7,8:14);
[mode,S] = get_mode_S(iso_results.x(1:7)-iso_results.x(8:14), sqrt(diag(Sp)));
nexttile; hold on;
plot(Ts, exp(iso_samples(:,:,1)-iso_samples(:,:,2)), Color=color_samp);
plot(Ts, polyval(true_params.gold_film.k,Ts)./polyval(true_params.gold_film.c,Ts)./polyval(true_params.gold_film.rho,Ts), '--', Color=color_true, LineWidth=lw);
errorbar(T, mode, mode-mode./S, mode.*S-mode, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);
xlim([290,910])
ylim([12.5,45])
xticks(300:200:900)
title("Gold-YSZ Experiment", Interpreter="latex")
ylabel("Film Diffusivity [mm$^2$/s]", Interpreter="latex")



indx = false(size(axisym_results.HessM,1),1);
indx(1:14) = true;
Hp = schur_comp(axisym_results.HessM, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');
Sp = Sp(1:7,1:7) + Sp(8:14,8:14) - 2*Sp(1:7,8:14);
sigma = sqrt(diag(Sp));
[mode, S] = get_mode_S(axisym_results.x(1:7)-axisym_results.x(8:14), sigma);
nexttile; hold on;
plot(Ts, exp(axisym_samples(:,:,1)-axisym_samples(:,:,2)), Color=color_samp);
plot(Ts, polyval(true_params.gold_film.k,Ts)./polyval(true_params.gold_film.c,Ts)./polyval(true_params.gold_film.rho,Ts), '--', Color=color_true, LineWidth=lw);
errorbar(T, mode, mode-mode./S, mode.*S-mode, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);
xlim([290,910])
ylim([12.5,45])
xticks(300:200:900)
title("Gold-Graphite Experiment", Interpreter="latex")



indx = false(size(iso_results.hessian,1),1);
indx(16:29) = true;
Hp = schur_comp(iso_results.hessian, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');
Sp = Sp(1:7,1:7) + Sp(8:14,8:14) - 2*Sp(1:7,8:14);
[mode,S] = get_mode_S(iso_results.x(16:22)-iso_results.x(23:29), sqrt(diag(Sp)));
nexttile; hold on;
plot(Ts, exp(iso_samples(:,:,3)-iso_samples(:,:,4)), Color=color_samp);
plot(Ts, polyval(true_params.ZrO2.k,Ts)./polyval(true_params.ZrO2.c,Ts)./polyval(true_params.ZrO2.rho,Ts), '--', Color=color_true, LineWidth=lw);
errorbar(T, mode, mode-mode./S, mode.*S-mode, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);
xlim([290,910])
xticks(300:200:900)
xlabel("Temperature [K]", Interpreter="latex")
ylabel("Substrate Diffusivity [mm$^2$/s]", Interpreter="latex")



nexttile; hold on;
plt_samp = plot(Ts, exp(axisym_samples(:,:,3)-axisym_samples(:,:,5)), Color=color_samp);
plot(Ts, exp(axisym_samples(:,:,4)-axisym_samples(:,:,5)), Color=color_samp);
plt_true = plot(Ts, polyval(true_params.graphite.k_perp,Ts)./polyval(true_params.graphite.c,Ts)./polyval(true_params.graphite.rho,Ts), '--', Color=color_true, LineWidth=lw);
plot(Ts, polyval(true_params.graphite.k_par,Ts)./polyval(true_params.graphite.c,Ts)./polyval(true_params.graphite.rho,Ts), '--', Color=color_true, LineWidth=lw);

indx = false(size(axisym_results.HessM,1),1);
indx([16:22,30:36]) = true;
Hp = schur_comp(axisym_results.HessM, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');
Sp = Sp(1:7,1:7) + Sp(8:14,8:14) - 2*Sp(1:7,8:14);
sigma = sqrt(diag(Sp));
[mode, S] = get_mode_S(axisym_results.x(16:22)-axisym_results.x(30:36), sigma);
plt_post = errorbar(T, mode, mode-mode./S, mode.*S-mode, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);

indx = false(size(axisym_results.HessM,1),1);
indx(23:36) = true;
Hp = schur_comp(axisym_results.HessM, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');
Sp = Sp(1:7,1:7) + Sp(8:14,8:14) - 2*Sp(1:7,8:14);
sigma = sqrt(diag(Sp));
[mode, S] = get_mode_S(axisym_results.x(23:29)-axisym_results.x(30:36), sigma);
errorbar(T, mode, mode-mode./S, mode.*S-mode, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);

text(560, 27, "$D_s^\perp$", Interpreter="latex")
plot([560,530],[27,27-3.55]-1, 'k')
text(460, 17.5, "$D_s^\parallel$", Interpreter="latex")
plot([460,430],[17.5,17.5-3.55]-1, 'k')

xlim([290,910])
xticks(300:200:900)
xlabel("Temperature [K]", Interpreter="latex")

lgd = legend([plt_post, plt_samp(1), plt_true], ["Posterior Mode and 95\% HDI$\quad\quad$", "Laplace Posterior Predictive Sample$\quad\quad$", "True Function Values$\quad\quad$"], Interpreter="latex", Orientation="horizontal");
lgd.Layout.Tile = 'north';

%% Correlation Matrices
bw = 0.05;
sc = 0.7;
gw = 0.25;
blw = 0.75;
grid_color = 0.5;
grid_bold = 0.1;
fig = figure;
fig.Units = 'inch';
fig.Position = [8 1 7 6.45];
addpath("drawbrace\")
tiledlayout(2,2,'TileSpacing','compact', "Padding","tight");

nexttile(1);
indx = false(size(iso_results.hessian,1),1);
indx([1:14,16:29]) = true;
Sp = inv(schur_comp(iso_results.hessian,indx));
Sp = 0.5*(Sp+Sp.');
imagesc(corrcov(Sp))
clim([-1,1])

D = size(Sp,1);
labels = "$\ln " + ["k_f", "C_f", "k_s", "C_s"]+"(\mathbf{T})$";
labels = labels(:);
yticks(4:7:D)
yticklabels(labels)
xticks(4:7:D)
xticklabels(labels)

ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0 0];
axis equal
xlim([0.5,D+0.5])
ylim([0.5,D+0.5])
title("Gold-YSZ Experiment", Interpreter="latex", FontSize=12)
ytickangle(90)
axis xy
ax.Clipping = 'off';
for i = 1:length(labels)
    drawbrace([D+0.5-7*(i-1),0.5],[D+0.5-7*i,0.5], bw, LineWidth=blw, Color='k');
    drawbrace([0.5, D+0.5-7*i],[0.5,D+0.5-7*(i-1)], bw, LineWidth=blw, Color='k');
end
ax.XAxis.TickLabelGapMultiplier = sc;
ax.YAxis.TickLabelGapMultiplier = sc;
hold on
plot((0.5:D).*ones(2,1), [0.5;D+0.5].*ones(1,D), Color=grid_color*ones(3,1), LineWidth=gw);
plot([0.5;D+0.5].*ones(1,D), (0.5:D).*ones(2,1), Color=grid_color*ones(3,1), LineWidth=gw);
plot((0.5:7:D).*ones(2,1), [0.5;D+0.5].*ones(1,D/7), Color=grid_bold*ones(3,1), LineWidth=gw);
plot([0.5;D+0.5].*ones(1,D/7), (0.5:7:D).*ones(2,1), Color=grid_bold*ones(3,1), LineWidth=gw);
nexttile(3);
A = [eye(7),-eye(7),zeros(7,7),zeros(7,7);
     zeros(7,7), zeros(7,7), eye(7), eye(7)];
S_lnD = A*Sp*A.';
D = size(S_lnD,1);
iso_corr = corrcov(S_lnD);
imagesc(iso_corr)
clim([-1,1])
c = colorbar;
addpath("BrewerMap-master\")
colormap(flipud(brewermap([],'RdBu')));
ylabel(c, "Correlation Coefficient", Interpreter="latex", FontSize=11)

labels = "$\ln " + ["D_f", "D_s"]+"(\mathbf{T})$";
labels = labels(:);
yticks(4:7:D)
yticklabels(labels)
xticks(4:7:D)
xticklabels(labels)
ytickangle(90)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0 0];
axis xy
ax.Clipping = 'off';
for i = 1:length(labels)
    drawbrace([D+0.5-7*(i-1),0.5],[D+0.5-7*i,0.5], bw, LineWidth=blw, Color='k');
    drawbrace([0.5, D+0.5-7*i],[0.5,D+0.5-7*(i-1)], bw, LineWidth=blw, Color='k');
end
ax.XAxis.TickLabelGapMultiplier = sc;
ax.YAxis.TickLabelGapMultiplier = sc;
hold on
plot((0.5:D).*ones(2,1), [0.5;D+0.5].*ones(1,D), Color=grid_color*ones(3,1), LineWidth=gw);
plot([0.5;D+0.5].*ones(1,D), (0.5:D).*ones(2,1), Color=grid_color*ones(3,1), LineWidth=gw);
plot((0.5:7:D).*ones(2,1), [0.5;D+0.5].*ones(1,D/7), Color=grid_bold*ones(3,1), LineWidth=gw);
plot([0.5;D+0.5].*ones(1,D/7), (0.5:7:D).*ones(2,1), Color=grid_bold*ones(3,1), LineWidth=gw);
axis equal
xlim([0.5,D+0.5])
ylim([0.5,D+0.5])
c.Layout.Tile = "east";

nexttile(2);

indx = false(size(axisym_results.HessM,1),1);
indx([1:14,16:36]) = true;
Sp = inv(schur_comp(axisym_results.HessM,indx));
Sp = 0.5*(Sp+Sp.');
D = size(Sp,1);
imagesc(corrcov(Sp))
clim([-1,1])

labels = "$\ln " + ["k_f", "C_f", "k_s^\perp", "k_s^\parallel", "C_s"]+"(\mathbf{T})$";
labels = labels(:);
yticks(4:7:D)
yticklabels(labels)
xticks(4:7:D)
xticklabels(labels)
ytickangle(90)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0 0];
title("Gold-Graphite Experiment", Interpreter="latex", FontSize=12)
axis xy
axis equal
xlim([0.5,D+0.5])
ylim([0.5,D+0.5])
xtickangle(0)
ax.Clipping = 'off';
for i = 1:length(labels)
    drawbrace([D+0.5-7*(i-1),0.5],[D+0.5-7*i,0.5], bw, LineWidth=blw, Color='k');
    drawbrace([0.5, D+0.5-7*i],[0.5,D+0.5-7*(i-1)], bw, LineWidth=blw, Color='k');
end
ax.XAxis.TickLabelGapMultiplier = sc;
ax.YAxis.TickLabelGapMultiplier = sc;
hold on
plot((0.5:D).*ones(2,1), [0.5;D+0.5].*ones(1,D), Color=grid_color*ones(3,1), LineWidth=gw);
plot([0.5;D+0.5].*ones(1,D), (0.5:D).*ones(2,1), Color=grid_color*ones(3,1), LineWidth=gw);
plot((0.5:7:D).*ones(2,1), [0.5;D+0.5].*ones(1,D/7), Color=grid_bold*ones(3,1), LineWidth=gw);
plot([0.5;D+0.5].*ones(1,D/7), (0.5:7:D).*ones(2,1), Color=grid_bold*ones(3,1), LineWidth=gw);
nexttile(4);
A = [ eye(7), -eye(7),  zeros(7),  zeros(7),  zeros(7);
      zeros(7), zeros(7),  eye(7), zeros(7), -eye(7);
      zeros(7), zeros(7), zeros(7),  eye(7), -eye(7)];
S_lnD = A*Sp*A.';
D = size(S_lnD,1);
axisym_corr = corrcov(S_lnD);
imagesc(axisym_corr)
clim([-1,1])
labels = "$\ln " + ["D_f", "D_s^\perp", "D_s^\parallel"]+"(\mathbf{T})$";
labels = labels(:);
yticks(4:7:D)
yticklabels(labels)
xticks(4:7:D)
xticklabels(labels)
ytickangle(90)
ax = gca;
ax.TickLabelInterpreter = 'latex';
axis equal
xlim([0.5,D+0.5])
ylim([0.5,D+0.5])
axis xy
ax.Clipping = 'off';
for i = 1:length(labels)
    drawbrace([D+0.5-7*(i-1),0.5],[D+0.5-7*i,0.5], bw, LineWidth=blw, Color='k');
    drawbrace([0.5, D+0.5-7*i],[0.5,D+0.5-7*(i-1)], bw, LineWidth=blw, Color='k');
end
ax.XAxis.TickLabelGapMultiplier = sc;
ax.YAxis.TickLabelGapMultiplier = sc;
hold on
plot((0.5:D).*ones(2,1), [0.5;D+0.5].*ones(1,D), Color=grid_color*ones(3,1), LineWidth=gw);
plot([0.5;D+0.5].*ones(1,D), (0.5:D).*ones(2,1), Color=grid_color*ones(3,1), LineWidth=gw);
plot((0.5:7:D).*ones(2,1), [0.5;D+0.5].*ones(1,D/7), Color=grid_bold*ones(3,1), LineWidth=gw);
plot([0.5;D+0.5].*ones(1,D/7), (0.5:7:D).*ones(2,1), Color=grid_bold*ones(3,1), LineWidth=gw);
%%
iso_corrs = diag(iso_corr(8:14,1:7));
axisym_corrs = [diag(axisym_corr(1:7,8:14)); diag(axisym_corr(1:7,15:21))];

iso_lnDf = iso_results.x(1:7)-iso_results.x(8:14);
iso_lnDs = iso_results.x(16:22)-iso_results.x(23:29);
axisym_lnDf = axisym_results.x(1:7)-axisym_results.x(8:14);
axisym_lnDs_perp = axisym_results.x(16:22)-axisym_results.x(30:36);
axisym_lnDs_par = axisym_results.x(23:29)-axisym_results.x(30:36);

figure
scatter(abs(iso_lnDf-iso_lnDs), iso_corrs, 50, "filled")
hold on;
scatter(abs(axisym_lnDf-axisym_lnDs_perp), axisym_corrs(1:7), "filled")
scatter(abs(axisym_lnDf-axisym_lnDs_par), axisym_corrs(8:14), "filled")
xlabel("Absolute Difference $\left|\ln D_f(T_i)-\ln D_s(T_i)\right|$", Interpreter="latex")
ylabel("Correlation Coefficient", Interpreter="latex")
legend( ...
    "Gold-YSZ", ...
    "Gold-Graphite (Perpendicular)$\;\;\;$", ...
    "Gold-Graphite (Parallel)", ...
    Interpreter="latex", Location="southeast", FontSize=10)
%% Refined GP Results (384x384)
iso_results = load("iso_MAP_results_refined_384.mat");

fprintf("ISO TABLE:\n\n")
for i = [15,30:34]
    indx = false(size(iso_results.hessian,1),1);
    indx(i) = true;
    Hp = schur_comp(iso_results.hessian, indx);
    mode = exp(iso_results.x(i)-1/Hp);
    S = get_logN_HDI_param(sigma=sqrt(1/Hp), alpha=0.95);
    fprintf("$\\mathrm{logN}(\\num{%.4g},\\num{%.4g})$ & \\num{%.4g} & [\\num{%.4g},\\num{%.4g}] \\\\\n", iso_results.x(i), 1/Hp, mode, mode/S, mode*S)
end
disp("\hline")
fprintf("\n\n")

% Iso Samples
indx = false(size(iso_results.hessian,1),1);
indx([1:14,16:29,31:34]) = true;
H_marg = schur_comp(iso_results.hessian, indx);
Sp = inv(H_marg);
Sp = 0.5*(Sp+Sp.');
mu = iso_results.x(indx);
mu0 = @(T) [iso_priors.priors.lnkf.mu, iso_priors.priors.lnCf.mu, iso_priors.priors.lnks.mu, iso_priors.priors.lnCs.mu] .* ones(length(T),1);
sigma0 = {...
    @(T) iso_priors.priors.lnkf.sigma .* ones(length(T),1);...
    @(T) iso_priors.priors.lnCf.sigma .* ones(length(T),1);...
    @(T) iso_priors.priors.lnks.sigma .* ones(length(T),1);...
    @(T) iso_priors.priors.lnCs.sigma .* ones(length(T),1);...
};
iso_samples = sample_posterior_predictive(T, Ts, mu0, sigma0, mu, Sp, N_samp);
iso_samples = reshape(iso_samples, N_samp, [], 4);

figure;
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

indx = false(size(iso_results.hessian,1),1);
indx(1:14) = true;
Hp = schur_comp(iso_results.hessian, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');
Sp = Sp(1:7,1:7) + Sp(8:14,8:14) - 2*Sp(1:7,8:14);
[mode,S] = get_mode_S(iso_results.x(1:7)-iso_results.x(8:14), sqrt(diag(Sp)));
nexttile; hold on;
plot(Ts, exp(iso_samples(:,:,1)-iso_samples(:,:,2)), Color=color_samp);
plot(Ts, polyval(true_params.gold_film.k,Ts)./polyval(true_params.gold_film.c,Ts)./polyval(true_params.gold_film.rho,Ts), '--', Color=color_true, LineWidth=lw);
errorbar(T, mode, mode-mode./S, mode.*S-mode, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);
xlim([290,910])
ylim([12.5,45])
xticks(300:200:900)
title("Gold-YSZ Experiment", Interpreter="latex")
ylabel("Film Diffusivity [mm$^2$/s]", Interpreter="latex")

indx = false(size(iso_results.hessian,1),1);
indx(16:29) = true;
Hp = schur_comp(iso_results.hessian, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');
Sp = Sp(1:7,1:7) + Sp(8:14,8:14) - 2*Sp(1:7,8:14);
[mode,S] = get_mode_S(iso_results.x(16:22)-iso_results.x(23:29), sqrt(diag(Sp)));
nexttile; hold on;
plt_samp = plot(Ts, exp(iso_samples(:,:,3)-iso_samples(:,:,4)), Color=color_samp);
plt_true = plot(Ts, polyval(true_params.ZrO2.k,Ts)./polyval(true_params.ZrO2.c,Ts)./polyval(true_params.ZrO2.rho,Ts), '--', Color=color_true, LineWidth=lw);
plt_post = errorbar(T, mode, mode-mode./S, mode.*S-mode, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);
xlim([290,910])
xticks(300:200:900)
xlabel("Temperature [K]", Interpreter="latex")
ylabel("Substrate Diffusivity [mm$^2$/s]", Interpreter="latex")

lgd = legend([plt_post, plt_samp(1), plt_true], ["Posterior Mode and 95\% HDI$\quad\quad$", "Laplace Posterior Predictive Sample$\quad\quad$", "True Function Values$\quad\quad$"], Interpreter="latex", Orientation="vertical");
lgd.Layout.Tile = 'north';

%%
iso_results = load("iso_MAP_results(true_kf_Cf).mat");
shift = 14;
% Iso Samples
indx = false(size(iso_results.hessian,1),1);
indx([16:29,33:34]-shift) = true;
H_marg = schur_comp(iso_results.hessian, indx);
Sp = inv(H_marg);
Sp = 0.5*(Sp+Sp.');
mu = iso_results.x(indx);

mu0 = @(T) [iso_priors.priors.lnks.mu, iso_priors.priors.lnCs.mu] .* ones(length(T),1);
sigma0 = {...
    @(T) iso_priors.priors.lnks.sigma .* ones(length(T),1);...
    @(T) iso_priors.priors.lnCs.sigma .* ones(length(T),1);...
};
iso_samples = sample_posterior_predictive(T, Ts, mu0, sigma0, mu, Sp, N_samp);
iso_samples = reshape(iso_samples, N_samp, [], 2);

indx = false(size(iso_results.hessian,1),1);
indx((16:29)-shift) = true;
Hp = schur_comp(iso_results.hessian, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');
Sp = Sp(1:7,1:7) + Sp(8:14,8:14) - 2*Sp(1:7,8:14);
[mode,S] = get_mode_S(iso_results.x((16:22)-shift)-iso_results.x((23:29)-shift), sqrt(diag(Sp)));
figure; hold on;
plt_samp = plot(Ts, exp(iso_samples(:,:,1)-iso_samples(:,:,2)), Color=color_samp);
plt_true = plot(Ts, polyval(true_params.ZrO2.k,Ts)./polyval(true_params.ZrO2.c,Ts)./polyval(true_params.ZrO2.rho,Ts), '--', Color=color_true, LineWidth=lw);
plt_post = errorbar(T, mode, mode-mode./S, mode.*S-mode, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);
xlim([290,910])
xticks(300:200:900)
title("Gold-YSZ Experiment", Interpreter="latex")
xlabel("Temperature [K]", Interpreter="latex")
ylabel("Substrate Diffusivity [mm$^2$/s]", Interpreter="latex")

lgd = legend([plt_post, plt_samp(1), plt_true], ["Posterior Mode and 95\% HDI$\quad\quad$", "Laplace Posterior Predictive Sample$\quad\quad$", "True Function Values$\quad\quad$"], Interpreter="latex", Orientation="vertical", Location="northoutside");%%
figure;
[X,Y,Z] = sphere(24);
X(Y<0) = NaN;
Z(Y<0) = NaN;
Y(Y<0) = NaN;
dx = 0.03;
Oi = [7, 4, 5, 8, 6];
load True_Params\gold_graphite_O.mat
Os_true = v(Oi,:);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for i = [3,5,1]
    nexttile;
    surf(X,zeros(size(X)),Z,FaceColor='w', EdgeColor=color_samp); hold on;
    
    indx = false(size(axisym_results.HessM,1),1);
    indx([36+i,41+i]) = true;
    H = schur_comp(axisym_results.HessM, indx);
    B = axisym_results.B([36+i,41+i,46+i],[36+i,41+i]);
    mu = axisym_results.x([36+i,41+i,46+i]);
    mu = mu/norm(mu);
    N = 361;
    pts2D = mvn95ellipse(H, zeros(1,2), N);
    pts3D = zeros(N,3);
    mu(3)
    for j = 1:N
        pts3D(j,:) = sphere_exp(mu, B*pts2D(:,j));
    end
    plot3(pts3D(:,1), zeros(N,1), pts3D(:,3), LineWidth=lw, Color=colors(i,:));
    plot3(pts3D(:,1), zeros(N,1), -pts3D(:,3), LineWidth=lw, Color=colors(i,:))
    scatter3(Os_true(i,1), 0, Os_true(i,3), 20, colors(i,:), LineWidth=lw, MarkerFaceColor='w', Marker="diamond");

    plot3(Os_true(i,1)+[-1,1,1,-1,-1]*dx, zeros(5,1), Os_true(i,3)+[-1,-1,1,1,-1]*dx, LineWidth=lw, Color='k');
    text(Os_true(i,1)-dx*1.2, 0, Os_true(i,3)-dx*1.1, "("+roman(i)+")", Interpreter="latex")

    axis off
    grid off
    axis equal
    view(0,180)

    xlim(Os_true(i,1)+[-1,1]*dx)
    zlim(Os_true(i,3)+[-1,1]*dx)
end

nexttile([2,2]);
plot(X, Z, 'Color', color_samp); hold on
plot(X.', Z.', 'Color', color_samp);
for i = 1:5
    indx = false(size(axisym_results.HessM,1),1);
    indx([36+i,41+i]) = true;
    H = schur_comp(axisym_results.HessM, indx);
    B = axisym_results.B([36+i,41+i,46+i],[36+i,41+i]);
    mu = axisym_results.x([36+i,41+i,46+i]);
    mu = mu/norm(mu);
    N = 361;
    pts2D = mvn95ellipse(H, zeros(1,2), N);
    pts3D = zeros(N,3);
    mu(3)
    for j = 1:N
        pts3D(j,:) = sphere_exp(mu, B*pts2D(:,j));
    end
    plot(pts3D(:,1), pts3D(:,3), LineWidth=lw, Color=colors(i,:));
    plot(pts3D(:,1), -pts3D(:,3), LineWidth=lw, Color=colors(i,:))
    scatter(Os_true(i,1), Os_true(i,3), 20, colors(i,:), LineWidth=lw, MarkerFaceColor='w', Marker="diamond");

    plot(Os_true(i,1)+[-1,1,1,-1,-1]*dx, Os_true(i,3)+[-1,-1,1,1,-1]*dx, LineWidth=lw, Color='k');
end

text(Os_true(1,1)-dx*3.6, Os_true(1,3)-dx*1.6, "("+roman(1)+")", Interpreter="latex")
text(Os_true(2,1)-dx*4, Os_true(2,3)-dx*1.6, "("+roman(2)+")", Interpreter="latex")
text(Os_true(3,1)-dx*4, Os_true(3,3)-dx*1.6, "("+roman(3)+")", Interpreter="latex")
text(Os_true(4,1)-dx*4.4, Os_true(4,3)-dx*1.6, "("+roman(4)+")", Interpreter="latex")
text(Os_true(5,1)-dx*4, Os_true(5,3)-dx*0.5, "("+roman(5)+")", Interpreter="latex")

xlabel("$x$", Interpreter="latex");
ylabel("$z$", Interpreter="latex")
grid off
axis equal
set(gca,'YDir','reverse')

[~,plt_HDI] = contour(zeros(2,2),zeros(2,2),zeros(2,2), LineWidth=lw, EdgeColor='k')
plt_true = scatter(2,2, 20, 'k', LineWidth=lw, MarkerFaceColor='w', Marker="diamond")
lgd = legend([plt_HDI, plt_true], ["Laplace Posterior 95\% HDI$\quad\quad$", "True Orientation$\quad\quad$"], Interpreter="latex", Orientation="horizontal", FontSize=10);
lgd.Layout.Tile = 'north';

xlim([-1,1])
zlim([-1,1])

for i = [2,4]
    nexttile;
    surf(X,zeros(size(X)),Z,FaceColor='w', EdgeColor=color_samp); hold on;
    
    indx = false(size(axisym_results.HessM,1),1);
    indx([36+i,41+i]) = true;
    H = schur_comp(axisym_results.HessM, indx);
    B = axisym_results.B([36+i,41+i,46+i],[36+i,41+i]);
    mu = axisym_results.x([36+i,41+i,46+i]);
    mu = mu/norm(mu);
    N = 361;
    pts2D = mvn95ellipse(H, zeros(1,2), N);
    pts3D = zeros(N,3);
    mu(3)
    for j = 1:N
        pts3D(j,:) = sphere_exp(mu, B*pts2D(:,j));
    end
    plot3(pts3D(:,1), zeros(N,1), pts3D(:,3), LineWidth=lw, Color=colors(i,:));
    plot3(pts3D(:,1), zeros(N,1), -pts3D(:,3), LineWidth=lw, Color=colors(i,:))
    scatter3(Os_true(i,1), 0, Os_true(i,3), 20, colors(i,:), LineWidth=lw, MarkerFaceColor='w', Marker="diamond");

    plot3(Os_true(i,1)+[-1,1,1,-1,-1]*dx, zeros(5,1), Os_true(i,3)+[-1,-1,1,1,-1]*dx, LineWidth=lw, Color='k');
    text(Os_true(i,1)-dx*1.2, 0, Os_true(i,3)-dx*1.1, "("+roman(i)+")", Interpreter="latex")

    axis off
    grid off
    axis equal
    view(0,180)

    xlim(Os_true(i,1)+[-1,1]*dx)
    zlim(Os_true(i,3)+[-1,1]*dx)
end

function [mode, S] = get_mode_S(mu, sigma)
    mode = exp(mu - sigma.^2);
    S = zeros(size(sigma));
    for i = 1:length(sigma)
        S(i) = get_logN_HDI_param(sigma=sigma(i), alpha=0.95);
    end
end

function str = roman(i)
    switch i
        case 1
            str = "c";
        case 2
            str = "d";
        case 3
            str = "a";
        case 4
            str = "e";
        case 5
            str = "b";
    end
end

function y = sphere_exp(x, v)
%SPHERE_EXP Exponential map on the unit n-sphere.
%
% y = sphere_exp(x,v)
%
% Inputs:
%   x : n-by-1 point on the unit sphere (norm(x)=1)
%   v : n-by-1 tangent vector at x (x'*v=0)
%
% Output:
%   y : n-by-1 point on the unit sphere

    % Ensure column vectors
    x = x(:);
    v = v(:);

    % Optional checks
    tol = 1e-10;
    if abs(norm(x) - 1) > tol
        error('x must have unit norm.');
    end
    if abs(x' * v) > tol
        error('v must lie in the tangent space: x''*v = 0.');
    end

    theta = norm(v);

    if theta < 1e-12
        y = x;
    else
        y = cos(theta)*x + sin(theta)/theta*v;
    end

    % Remove roundoff error
    y = y / norm(y);

end

function X = mvn95ellipse(Lambda, mu, N)
% Generate N points on the 95% confidence ellipse
% of a 2D Gaussian specified by its precision matrix.
%
% Lambda : 2x2 precision matrix
% mu     : 2x1 mean vector
% N      : number of points (default 100)

    if nargin < 3
        N = 100;
    end

    mu = mu(:);

    % 95% chi-square quantile for 2 DOF
    c = 5.991464547;

    % Cholesky factor of precision matrix
    L = chol(Lambda,'lower');

    % Unit circle
    theta = linspace(0,2*pi,N);
    U = [cos(theta); sin(theta)];

    % Transform
    X = mu + sqrt(c) * (L'\U);

end