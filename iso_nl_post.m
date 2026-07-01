function [psi, grad] = iso_nl_post(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, tau, theta, kD, obs, lnsx, lnsy, f, Nx, dr, T, priors)
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
        psi_theta    = nl_theta_prior(theta, priors.theta.mu, priors.theta.sigma);
        phi_GP_lnkf  = nl_GP_prior(lnkf, priors.lnkf.mu*ones(size(lnkf)), priors.lnkf.sigma*ones(size(lnkf)), T, theta(1));
        phi_GP_lnCf  = nl_GP_prior(lnCf, priors.lnCf.mu*ones(size(lnCf)), priors.lnCf.sigma*ones(size(lnCf)), T, theta(2));
        phi_GP_lnhf  = nln(lnhf, priors.lnhf.mu, priors.lnhf.sigma, [true, false, false]);
        phi_GP_lnks  = nl_GP_prior(lnks, priors.lnks.mu*ones(size(lnks)), priors.lnks.sigma*ones(size(lnks)), T, theta(3));
        phi_GP_lnCs  = nl_GP_prior(lnCs, priors.lnCs.mu*ones(size(lnCs)), priors.lnCs.sigma*ones(size(lnCs)), T, theta(4));
        psi_tau = nln(tau, priors.tau.mu, priors.tau.sigma, [true,false,false]);
    else
        len_M_vars = length(lnkf) + length(lnCf) + length(lnaf) + length(lnhf) + length(lnks) + length(lnCs) + length(lnas) + length(lnRth);
        len_kT = length(kT);
        len_theta = length(theta);
    
        [psi_nll, grad_nll] = iso_nll_M(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, kT, kD, obs, lnsx, lnsy, f, Nx, dr);
        grad_nll = [grad_nll; zeros(len_theta, 1)];
    
        [psi_theta, grad_theta] = nl_theta_prior(theta, priors.theta.mu, priors.theta.sigma);
        grad_theta = [zeros(len_M_vars + len_kT,1); grad_theta];
    
        [phi_GP_lnkf, grad_GP_lnkf, grad_GP_theta1] = nl_GP_prior(lnkf, priors.lnkf.mu*ones(size(lnkf)), priors.lnkf.sigma*ones(size(lnkf)), T, theta(1));
        [phi_GP_lnCf, grad_GP_lnCf, grad_GP_theta2] = nl_GP_prior(lnCf, priors.lnCf.mu*ones(size(lnCf)), priors.lnCf.sigma*ones(size(lnCf)), T, theta(2));
        [phi_GP_lnhf, grad_GP_lnhf]   = nln(lnhf, priors.lnhf.mu, priors.lnhf.sigma, [true, false, false]);
        [phi_GP_lnks, grad_GP_lnks, grad_GP_theta3] = nl_GP_prior(lnks, priors.lnks.mu*ones(size(lnks)), priors.lnks.sigma*ones(size(lnks)), T, theta(3));
        [phi_GP_lnCs, grad_GP_lnCs, grad_GP_theta4] = nl_GP_prior(lnCs, priors.lnCs.mu*ones(size(lnCs)), priors.lnCs.sigma*ones(size(lnCs)), T, theta(4));
        [psi_tau, grad_tau] = nln(tau, priors.tau.mu, priors.tau.sigma, [true,false,false]);

        grad_GP = [grad_GP_lnkf; grad_GP_lnCf; 0; grad_GP_lnhf{1}; grad_GP_lnks; grad_GP_lnCs; 0; 0; 0; grad_GP_theta1; grad_GP_theta2; grad_GP_theta3; grad_GP_theta4];

        grad = grad_nll + grad_theta + grad_GP;
        grad(33) = kT*grad(33) + grad_tau{1};
        grad = grad([1:14,16:30,33:37]);
        % if any(isnan(grad))
        %     pause(1);
        % end
    end
    psi = psi_nll + psi_theta + phi_GP_lnkf + phi_GP_lnCf + phi_GP_lnhf + phi_GP_lnks + phi_GP_lnCs + psi_tau;
    % if isnan(psi)
    %     pause(1);
    % end
end