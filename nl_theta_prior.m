function [psi, grad] = nl_theta_prior(theta, mu_theta, sigma_theta)
    [psi, grad] = nln(theta,mu_theta,sigma_theta,[true, false, false]);
    psi = sum(psi);
    grad = grad{1};
end