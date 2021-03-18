function [x, n, error] = GMRES(A, b, maximum_iteration_number, tol)
    %initial arbitrary vector q_0
    q_1 = b / norm(b);
    n = length(A);
    H_tilda = zeros(n+1, n);
    Q = zeros(length(q_1), n+1);
    e_1 = zeros(n+1, 1);
    e_1(1,1) = 1;
    x = 0;
    Q(:,1) = q_1;
    error = zeros(1, maximum_iteration_number);
    for n = 1:maximum_iteration_number
       %step n of arnoldi iteration
       %v = A*q_n;
       v = A*Q(:,n);
       for j=1:n
          %h_j_n = q_j' * v;
          H_tilda(j,n) = Q(:,j)' * v;
          %v = v - h_j_n * q_j;
          v = v - H_tilda(j,n) * Q(:,j);
       end
       %h_n_next_n = norm(v);
       H_tilda(n+1,n) = norm(v);
       %q_n_next = v / h_n_next_n;
       Q(:,n+1) = v / H_tilda(n+1,n);
       %end of step n of Arnoldi iteration
       
       %find y to minimize ||r_n|| using QR
       [Q_r, R] = qr(H_tilda(1:n+1, 1:n));
       d = Q_r' * (norm(b)*e_1(1:n+1));
       y = R\d;
       x = Q(:, 1:n)*y;
       error(n) = norm(A*x - b);
       if error(n) < tol
           return
       end
    end
end

