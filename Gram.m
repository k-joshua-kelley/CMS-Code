function [K, grad_l] = Gram(X,Xp,sigma,l)
arguments
    X     (:,1) double
    Xp    (:,1) double
    sigma (:,1) double
    l     (1,1) double
end

% pairwise squared distances
r2 = pdist2(X, Xp, 'squaredeuclidean');

% base kernel
E = exp(-r2/(2*l^2));
K = (sigma * sigma.') .* E;

if nargout > 1
    % dK/dl
    grad_l = K .* (r2 / l^3);
end
end