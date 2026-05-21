function [psi, grad] = nl_Bigham_prior_3D(O, Q, k)
    psi = 0;
    grad = zeros(size(O));
    for i = 1:size(O,1)
        [a, b] = unnorm_nl_Bingham(O(i,:), Q, k);
        psi = psi + a;
        grad(i,:) = b;
    end
    grad = grad(:);
end