function [psi, grad] = unnorm_nl_Bingham(x, Q, k)
arguments
    x (:,1) double
    Q (:,:) double
    k (:,1) double
end

    A = Q * diag([0; k]) * Q.';
    
    psi = x.' * A * x;

    if nargout > 1
        grad = 2 * A * x;
    end

end