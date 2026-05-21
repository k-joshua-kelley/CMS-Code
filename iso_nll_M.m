function [pnll, g_nll_mv_lnkT] = iso_nll_M(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, kT, kD, obs, lnsx, lnsy, f, Nx, dr)
    obs = permute(obs, [4,1,2,3]);
    kD = permute(kD, [4,1,2,3]);

    if nargout < 2
        phi_pred = iso_fm_phi(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, lnsx, lnsy, f, Nx, dr);
        pnll = nll(phi_pred, obs, kT, kD);
    else
        [phi_pred, J_fm] = iso_fm_phi(lnkf, lnCf, lnaf, lnhf, lnks, lnCs, lnas, lnRth, lnsx, lnsy, f, Nx, dr);

        [pnll, g_nll_phi, g_nll_kT] = nll(phi_pred, obs, kT, kD);
        g_nll_mv  = J_fm.' * g_nll_phi(:);
        g_nll_kT = sum(g_nll_kT, "all");
        
        g_nll_mv_lnkT = [g_nll_mv; g_nll_kT];
    end

    pnll = sum(pnll, "all");
end