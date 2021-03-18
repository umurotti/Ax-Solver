function [S] = produce_S(m,tau)
    %put 1 at each diagonal position
    d = ones(m,1); % The diagonal values
    %random number from a unifrom distrubtion on [-1, 1]
    a = -1;
    b = 1;
    t = triu(bsxfun(@min,d,d.').*(a + (b-a).*rand(m)),1); % The upper trianglar random values
    S = diag(d)+t+t.'; % Put them together in a symmetric matrix
    %
    S(abs(S) > tau) = 0;
    S(find(eye(m))) = 1;
end

