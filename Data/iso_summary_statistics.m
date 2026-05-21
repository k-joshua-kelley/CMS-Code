clear; clc; close all;

load iso_noisy_data.mat

[N, n, Nf, Nprobe, Nrep] = size(iso_phi_noisy);
options = optimoptions("fmincon");
mu = zeros(N,n,Nf,Nprobe);
kappa = zeros(N,n,Nf,Nprobe);
for i = 1:N
    for j = 1:n
        for k = 1:Nf
            for l = 1:Nprobe
                x = fmincon(@(x) sum(vM(iso_phi_noisy(i,j,k,l,:), x(1), x(2))), [0, 100], [],[],[],[],[-pi,0],[pi,inf]);
                mu(i,j,k,l) = x(1);
                kappa(i,j,k,l) = x(2);
            end
        end
    end
end

for i = randperm(N,3)
    j = 1;
    for k = randperm(Nf,3)
        for l = randperm(Nprobe,3)
            clf
            histogram(iso_phi_noisy(i,j,k,l,:), Normalization="pdf")
            hold on;
            x = linspace(min(iso_phi_noisy(i,j,k,l,:)), max(iso_phi_noisy(i,j,k,l,:)), 201);
            plot(x, exp(-vM(x,mu(i,j,k,l), kappa(i,j,k,l))), LineWidth=2)
            drawnow;
            pause(1)
        end
    end
end

save("iso_data_stats.mat", "mu", "kappa", "T", "f", "dr")