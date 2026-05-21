clear; clc; close all;

load axisym_noisy_data.mat

[N, n, Nf, Nprobe, Nrep] = size(axisym_phi_noisy);
options = optimoptions("fmincon");
mu = zeros(N,n,Nf,Nprobe);
kappa = zeros(N,n,Nf,Nprobe);
count = 0;
for i = 1:N
    for j = 1:n
        for k = 1:Nf
            for l = 1:Nprobe
                [x, ~, exitflag] = fmincon(@(x) sum(vM(axisym_phi_noisy(i,j,k,l,:), x(1), x(2))), [0, 300], [],[],[],[],[-pi,0],[pi,inf]);
                if exitflag ~= 1 && exitflag ~= 2
                    pause(1);
                end
                mu(i,j,k,l) = x(1);
                kappa(i,j,k,l) = x(2);
                count = count + 1;
                disp("Progress: " + count + "/" + (N*n*Nf*Nprobe))
            end
        end
    end
end

save("axisym_data_stats.mat", "mu", "kappa", "T", "f", "X_probe")