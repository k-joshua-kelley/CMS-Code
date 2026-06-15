function H = central_diff_hessian(fun, x, h)
    n = length(x);
    H = zeros(n,n);
    for j = 1:n
        x0 = x;
        x0(j) = x0(j) + h;
        [~,gp] = fun(x0);
        x0 = x;
        x0(j) = x0(j) - h;
        [~,gm] = fun(x0);
        H(:,j) = 0.5*(gp+gm)/h;
        assert(length(gp)==n);
        assert(length(gm)==n);
    end
    H = 0.5*(H+H.');
end