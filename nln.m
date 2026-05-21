function [psi,grad] = nln(x,mu,sigma,grad_indx)
arguments
    x (:,1) double
    mu (:,1) double
    sigma (:,1) double
    grad_indx (1,3) logical = true(1,3);
end

% --- size check ---
lens = [numel(x), numel(mu), numel(sigma)];
lens = lens(lens > 1);
assert(isempty(lens) || all(lens == lens(1)), ...
    'Non-scalar inputs must have the same length.');

% --- core computation ---
diff = x - mu;
inv_sigma = 1 ./ sigma;
inv_sigma2 = inv_sigma.^2;

psi = 0.5 * (diff.^2 .* inv_sigma2) + log(sigma) + 0.5*log(2*pi);

% --- gradients ---
if nargout > 1
    grad = cell(1,3);

    if grad_indx(1)  % d/dx
        grad{1} = diff .* inv_sigma2;
    end

    if grad_indx(2)  % d/dmu
        grad{2} = -diff .* inv_sigma2;
    end

    if grad_indx(3)  % d/dsigma
        grad{3} = inv_sigma - (diff.^2) .* (inv_sigma.^3);
    end
end
end