function [psi, grad] = axisym_nl_post(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD, theta, T, priors)
    if nargout < 2
        psi_nll       = axisym_nll_M(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD);
        psi_theta     = nl_theta_prior(theta, priors.theta.mu, priors.theta.sigma);
        psi_lnkf      = nl_GP_prior(lnkf, priors.lnkf.mu*ones(size(lnkf)), priors.lnkf.sigma*ones(size(lnkf)), T, theta(1));
        psi_lnCf      = nl_GP_prior(lnCf, priors.lnCf.mu*ones(size(lnCf)), priors.lnCf.sigma*ones(size(lnCf)), T, theta(2));
        psi_lnhf      = nln(lnhf, priors.lnhf.mu, priors.lnhf.sigma, [true, false, false]);
        psi_lnks_perp = nl_GP_prior(lnks_perp, priors.lnks_perp.mu*ones(size(lnks_perp)), priors.lnks_perp.sigma*ones(size(lnks_perp)), T, theta(3));
        psi_lnks_par  = nl_GP_prior(lnks_par, priors.lnks_par.mu*ones(size(lnks_par)), priors.lnks_par.sigma*ones(size(lnks_par)), T, theta(4));
        psi_lnCs      = nl_GP_prior(lnCs, priors.lnCs.mu*ones(size(lnCs)), priors.lnCs.sigma*ones(size(lnCs)), T, theta(5));
        psi_Os        = nl_Bigham_prior_3D([Os1(:), Os2(:), Os3(:)], priors.Os.Q, priors.Os.kappas);
    else
        len_M_vars = length(lnkf) + length(lnCf) + length(lnhf) + length(lnks_perp) + length(lnks_par) + length(lnCs);
        len_O = length(Os1);
        len_kT = length(kT);
        len_theta = length(theta);
    
        [psi_nll, grad_nll] = axisym_nll_M(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD);
        grad_nll = [grad_nll; zeros(len_theta, 1)];
    
        [psi_theta, grad_theta] = nl_theta_prior(theta, priors.theta.mu, priors.theta.sigma);
        grad_theta = [zeros(len_M_vars + len_kT + 3*len_O,1); grad_theta];

        [psi_lnkf, grad_GP_lnkf, grad_GP_theta1] = nl_GP_prior(lnkf, priors.lnkf.mu*ones(size(lnkf)), priors.lnkf.sigma*ones(size(lnkf)), T, theta(1));
        [psi_lnCf, grad_GP_lnCf, grad_GP_theta2] = nl_GP_prior(lnCf, priors.lnCf.mu*ones(size(lnCf)), priors.lnCf.sigma*ones(size(lnCf)), T, theta(2));
        [psi_lnhf, grad_GP_lnhf] = nln(lnhf, priors.lnhf.mu, priors.lnhf.sigma, [true, false, false]);
        [psi_lnks_perp, grad_GP_lnks_perp, grad_GP_theta3] = nl_GP_prior(lnks_perp, priors.lnks_perp.mu*ones(size(lnks_perp)), priors.lnks_perp.sigma*ones(size(lnks_perp)), T, theta(3));
        [psi_lnks_par, grad_GP_lnks_par, grad_GP_theta4] = nl_GP_prior(lnks_par, priors.lnks_par.mu*ones(size(lnks_par)), priors.lnks_par.sigma*ones(size(lnks_par)), T, theta(4));
        [psi_lnCs, grad_GP_lnCs, grad_GP_theta5] = nl_GP_prior(lnCs, priors.lnCs.mu*ones(size(lnCs)), priors.lnCs.sigma*ones(size(lnCs)), T, theta(5));
        [psi_Os, grad_Os] = nl_Bigham_prior_3D([Os1(:), Os2(:), Os3(:)], priors.Os.Q, priors.Os.kappas);

        grad_GP = [grad_GP_lnkf; grad_GP_lnCf; grad_GP_lnhf{1}; grad_GP_lnks_perp; grad_GP_lnks_par; grad_GP_lnCs; grad_Os; 0; grad_GP_theta1; grad_GP_theta2; grad_GP_theta3; grad_GP_theta4; grad_GP_theta5];

        grad = grad_nll + grad_theta + grad_GP;
    end
    psi = psi_nll + psi_theta + psi_lnkf + psi_lnCf + psi_lnhf + psi_lnks_perp + psi_lnks_par + psi_lnCs + psi_Os;
end