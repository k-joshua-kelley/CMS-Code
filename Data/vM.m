function [p, g, H] = vM(x, m, k, vars)
arguments
    x (:,1) double {mustBeReal}
    m (:,1) double {mustBeReal}
    k (:,1) double {mustBeNonnegative}
    vars (:,1) = [true; true; true]
end
    c = cos(x-m);
    kc = k.*c;
    I0 = besseli(0,k,1);
    p = -kc + log(2*pi*I0) + k;
    if nargout > 1
        s = sin(x-m);
        ks = k.*s;
        g = zeros(length(p), 3);
        mask = false(3,1);
        mask(vars) = true;
        if any(mask(1:2))
            g(:,1) = -ks;
            g(:,2) = ks;
        end
        if mask(3)
            I1_I0 = besseli(1,k,1) ./ I0;
            g(:,3) = I1_I0 - c;
        end
        if nargout > 2
            H = zeros(length(p), 3, 3);
            if any(mask(1:2))
                H(:,1,1) = kc; H(:,1,2) = -kc; H(:,1,3) = -s;
                               H(:,2,2) =  kc; H(:,2,3) =  s;
            end
            if mask(3)
                H(:,3,3) = 1 - I1_I0./k - I1_I0.*I1_I0;
            end
        end
    end
    if nargout > 1
    g = g(mask);
    end
    if nargout > 2
    H = H(mask);
    end
end