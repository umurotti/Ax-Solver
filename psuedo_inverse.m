function [A_psuedo] = psuedo_inverse(A)
    [U,S,V] = svd(A);
    A_psuedo = V * transpose(diag(1 ./ (S*ones(size(S, 1),1)))) * transpose(U);
end

