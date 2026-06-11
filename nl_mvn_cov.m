function [psi,grad] = nl_mvn_cov(x,mu,C,grad_indx)
arguments (Input)
    x  (:,1) double
    mu (:,1) double
    C  (:,:) double
    grad_indx (1,3) logical = true(1,3)
end
arguments (Output)
    psi  (1,1) double
    grad (:,1) cell
end

d = length(x);
diff = x - mu;

% Cholesky
L = chol(C+1e-6*eye(size(C)),'lower');

% log det
logdetS = 2*sum(log(diag(L)));

% quadratic
y = L \ diff;
quad = sum(y.^2);

% negative log pdf
psi = 0.5*d*log(2*pi) + 0.5*logdetS + 0.5*quad;

if nargout > 1

    % S^{-1} * diff
    Sinv_diff = L' \ y;

    % Precompute S^{-1}
    Sinv = L' \ (L \ eye(d));

    grad = {};

    % grad w.r.t x
    if grad_indx(1)
        grad{end+1} = Sinv_diff;
    end

    % grad w.r.t mu
    if grad_indx(2)
        grad{end+1} = -Sinv_diff;
    end

    % grad w.r.t C (match parameterization!)
    if grad_indx(3)
        grad{end+1} = 0.5 * (Sinv - Sinv*(diff*diff.')*Sinv);
    end
end

end
