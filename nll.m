function [p,g_x, g_kT] = nll(x, mu, kT, kD)
    diff = x - mu;
    c = cos(diff);
    kT2 = kT.*kT; kD2 = kD.*kD; kTkD = kT.*kD;
    kS = sqrt(kT2 + kD2 + 2*kTkD.*c);
    I0Stilde = besseli(0,kS,1);
    log_I0S = log(I0Stilde) + kS;

    norm_const = log(besseli(0,kT,1)) + kT + log(besseli(0,kD,1)) + kD + log(2*pi);

    p = norm_const - log_I0S;

    if nargout > 1
        s = sin(diff);
        kTkD_kS = kTkD./kS;
        I1I0S = besseli(1,kS,1) ./ I0Stilde;

        g_x = I1I0S .* kTkD_kS .* s;
    end

    if nargout > 2
        % ratio I1/I0 at kT (use scaled for stability)
        I1I0T = besseli(1,kT,1) ./ besseli(0,kT,1);
    
        % reuse existing terms
        dkS_dkT = (kT + kD.*c) ./ kS;
    
        g_kT = I1I0T - I1I0S .* dkS_dkT;
    end
end