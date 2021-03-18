function [x, n, error] = MINRES(A, b, maximum_iteration_number, tol)
    %initial arbitrary vector q_0
    q_0 = 0;
    q_n = q_0;
    q_1 = b / norm(b);
    beta_0 = 0;
    beta_n = beta_0;
    n = length(A);
    Q = zeros(length(q_1), n+1);
    T_tilda = zeros(n+1, n);
    e_1 = zeros(n+1, 1);
    e_1(1,1) = 1;
    x = 0;
    Q(:,1) = q_1;
    error = zeros(1, maximum_iteration_number);
    for n = 1:maximum_iteration_number
        beta_n_prev = beta_n;
        q_n_prev = q_n;
        %Lancsoz step
        %v = Aqn
        v = A*Q(:,n);
        alpha_n = Q(:,n)'*v;
        v = v - beta_n_prev*q_n_prev-alpha_n*Q(:,n);
        beta_n = norm(v);
        Q(:,n+1) = v/beta_n;
        T_tilda(n,n) = alpha_n;
        T_tilda(n+1,n) = beta_n;
        T_tilda(n,n+1) = beta_n;
        %end of step n of Lancsoz iteration
        
        %find y to minimize ||r_n|| using QR
        [Q_r, R] = qr(T_tilda(1:n+1, 1:n));
        d = Q_r' * (norm(b)*e_1(1:n+1));
        y = R\d;
        
        x = Q(:, 1:n)*y;
        error(n) = norm(A*x - b);
        if error(n) < tol
            return
        end
        q_n = Q(:,n);
    end
end

