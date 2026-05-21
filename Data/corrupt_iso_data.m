clear
clc
close all
rng(22)

load("COMSOL_iso_results/COMSOL_iso_results.mat")

dr = [5,10,20];
Nrep = 100;
bias_sigma = deg2rad(1);

[N,n,Nf,~] = size(T0tilde);
iso_phi = zeros(N,n,Nf,length(dr));

for i = 1:length(dr)
    [~,indx] = min(abs(dr(i)-r));
    iso_phi(:, :, :, i) = angle(T0tilde(:,:,:,indx));
end

T = (300:100:900).';

bias = normrnd(0,bias_sigma,size(iso_phi));
noise_sigma = deg2rad(rand(size(iso_phi))*2+1);

iso_phi_noisy = zeros([size(iso_phi),Nrep]);
for ri = 1:Nrep
    iso_phi_noisy(:,:,:,:,ri) = iso_phi + bias + normrnd(0,noise_sigma,size(iso_phi));
end

save("iso_noisy_data.mat","iso_phi_noisy", "T", "f", "dr")
save("iso_noise_free_data.mat","iso_phi", "T", "f", "dr")