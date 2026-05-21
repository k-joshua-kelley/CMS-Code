function [pnll, g_nll_mv_lnkT] = axisym_nll_M(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max, obs, kT, kD)
    if nargout < 2
        phi_pred = axisym_fm_phi(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max);
        pnll = nll(phi_pred, obs, kT, kD);
    else
        [phi_pred, J_fm] = axisym_fm_phi(lnkf,lnCf,lnaf,lnhf,lnks_perp,lnks_par,lnCs,lnas,lnRth,Os1,Os2,Os3,lnsx,lnsy,f,Nx,X_probe,x_max);

        [pnll, g_nll_phi, g_nll_kT] = nll(phi_pred, obs, kT, kD);
        if numel(g_nll_phi) ~= size(J_fm,1)
            g_nll_phi = sum(g_nll_phi,2);
        end
        g_nll_mv  = J_fm.' * g_nll_phi(:);
        g_nll_kT = sum(g_nll_kT, "all");
        
        g_nll_mv_lnkT = [g_nll_mv; g_nll_kT];
    end

    pnll = sum(pnll, "all");
end