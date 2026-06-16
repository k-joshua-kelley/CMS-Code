function H = central_diff_hessian(fun, x, h)
    n = length(x);
    H = zeros(n,n);

    for j = 1:n
        xp = x; xp(j) = xp(j) + h;
        xm = x; xm(j) = xm(j) - h;

        [~,gp] = fun(xp);
        [~,gm] = fun(xm);

        save(j + "_gradients.mat", "xp", "gp", "xm", "gm")

        assert(length(gp)==n);
        assert(length(gm)==n);

        H(:,j) = (gp - gm) / (2*h);

        disp("Progress:" + j + "/" + n)
    end
    save("raw_H.mat", "H")

    H = 0.5*(H + H.');
end