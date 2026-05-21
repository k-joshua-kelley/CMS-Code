clear; clc; close all; rng(22);

COMSOL_results = load("COMSOL_axisym_results/COMSOL_axisym_results.mat");

dr = [5,10,20];
sweep_angle = deg2rad(0:15:90).';
x_probe = cos(sweep_angle(:)).*dr(:).';
y_probe = sin(sweep_angle(:)).*dr(:).';
X_probe = [x_probe(:), y_probe(:)];

Nrep = 100;
bias_sigma = deg2rad(1);

[N,n,Nf,~] = size(COMSOL_results.T0tilde);
Nprobe = size(X_probe,1);
axisym_phi = zeros(N,n,Nf,Nprobe);

for i = 1:Nprobe
    [~,indx] = min(pdist2(X_probe(i,:), [COMSOL_results.x_probe(:), COMSOL_results.y_probe(:)], "euclidean"));
    axisym_phi(:, :, :, i) = angle(COMSOL_results.T0tilde(:,:,:,indx));
end

T = (300:100:900).';
f = COMSOL_results.f;


bias = normrnd(0,bias_sigma,size(axisym_phi));
noise_sigma = deg2rad(rand(size(axisym_phi))*2+1);

axisym_phi_noisy = zeros([size(axisym_phi),Nrep]);
for ri = 1:Nrep
    axisym_phi_noisy(:,:,:,:,ri) = axisym_phi + bias + normrnd(0,noise_sigma,size(axisym_phi));
end

save("axisym_noisy_data.mat","axisym_phi_noisy", "T", "f", "X_probe")
save("axisym_noise_free_data.mat","axisym_phi", "T", "f", "X_probe")