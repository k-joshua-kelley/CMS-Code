function [psi, grad_x, grad_theta] = nl_GP_prior(x, mu, sigma, X, theta)
arguments
    x     (:,1)
    mu    (:,1)
    sigma (:,1)
    X     (:,1)
    theta (1,1)
end
lens = [numel(x), numel(mu), numel(X), numel(X), numel(sigma)];
lens = lens(lens > 1);
assert(isempty(lens) || all(lens == lens(1)), ...
    'Non-scalar inputs must have the same length.');

l = exp(theta);
if nargout < 2
    psi = nl_mvn_cov(x, mu, Gram(X,X,sigma,l)+1e-6*eye(length(X)));
else
    [K, grad_K_l] = Gram(X,X,sigma,l);
    grad_K_theta = l * grad_K_l;
    
    [psi, grad_cell] = nl_mvn_cov(x, mu, K+1e-6*eye(length(X)), [true, false, true]);
    grad_x = grad_cell{1};
    grad_theta = sum(grad_cell{2} .* grad_K_theta, "all");
end
