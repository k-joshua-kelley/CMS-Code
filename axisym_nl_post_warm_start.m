function [psi, grad] = axisym_nl_post_warm_start(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD, priors)
arguments
    lnkf (1,1)
    lnCf (1,1) 
    lnaf (1,1)
    lnhf (1,1)
    lnks_perp (1,1)
    lnks_par (1,1)
    lnCs (1,1)
    lnas (1,1)
    lnRth (1,1)
    Os1 (1,1)
    Os2 (1,1)
    Os3 (1,1)
    lnsx (1,1)
    lnsy (1,1)
    f (:,1)
    Nx (1,1)
    X_probe (:,2)
    x_max
    obs (1,:,:,:)
    kT (1,1)
    kD (1,:,:,:)
    priors 
end
    if nargout < 2
        psi_nll       = axisym_nll_M(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD);
        psi_lnkf      = nln(lnkf, priors.lnkf.mu, priors.lnkf.sigma, [true, false, false]);
        psi_lnCf      = nln(lnCf, priors.lnCf.mu, priors.lnCf.sigma, [true, false, false]);
        psi_lnhf      = nln(lnhf, priors.lnhf.mu, priors.lnhf.sigma, [true, false, false]);
        psi_lnks_perp = nln(lnks_perp, priors.lnks_perp.mu, priors.lnks_perp.sigma, [true, false, false]);
        psi_lnks_par  = nln(lnks_par, priors.lnks_par.mu, priors.lnks_par.sigma, [true, false, false]);
        psi_lnCs      = nln(lnCs, priors.lnCs.mu, priors.lnCs.sigma, [true, false, false]);
        psi_Os        = nl_Bigham_prior_3D([Os1, Os2, Os3], priors.Os.Q, priors.Os.kappas);
    else
        [psi_nll, grad_nll] = axisym_nll_M(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD);
        grad_nll = grad_nll(1:end-1);

        [psi_lnkf,grad_lnkf]      = nln(lnkf, priors.lnkf.mu, priors.lnkf.sigma, [true, false, false]);
        [psi_lnCf,grad_lnCf]      = nln(lnCf, priors.lnCf.mu, priors.lnCf.sigma, [true, false, false]);
        [psi_lnhf,grad_lnhf]      = nln(lnhf, priors.lnhf.mu, priors.lnhf.sigma, [true, false, false]);
        [psi_lnks_perp,grad_lnks_perp] = nln(lnks_perp, priors.lnks_perp.mu, priors.lnks_perp.sigma, [true, false, false]);
        [psi_lnks_par,grad_lnks_par]  = nln(lnks_par, priors.lnks_par.mu, priors.lnks_par.sigma, [true, false, false]);
        [psi_lnCs,grad_lnCs]      = nln(lnCs, priors.lnCs.mu, priors.lnCs.sigma, [true, false, false]);
        [psi_Os,grad_Os]        = nl_Bigham_prior_3D([Os1, Os2, Os3], priors.Os.Q, priors.Os.kappas);

        grad_GP = [grad_lnkf{1}; grad_lnCf{1}; grad_lnhf{1}; grad_lnks_perp{1}; grad_lnks_par{1}; grad_lnCs{1}; grad_Os];

        grad = grad_nll + grad_GP;
    end
    psi = psi_nll + psi_lnkf + psi_lnCf + psi_lnhf + psi_lnks_perp + psi_lnks_par + psi_lnCs + psi_Os;
end