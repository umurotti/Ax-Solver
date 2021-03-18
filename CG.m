function [x_n, i, error] = CG(A, b, maximum_iteration_number, tol)
%     x_0 = 0;
%     r_0 = b;
%     p_0 = r_0;
%     
%     x_n = x_0;
%     r_n = r_0;
%     p_n = p_0;
    
    x_n = 0;
    r_n = b;
    p_n = r_n;
    error = zeros(1, maximum_iteration_number);
    for i = 1:maximum_iteration_number
        %update
        r_n_prev = r_n;
        p_n_prev = p_n;
        x_n_prev = x_n;
        %
        alpha_n = (transpose(r_n_prev)*r_n_prev)/(transpose(p_n_prev)*A*p_n_prev);
        x_n = x_n_prev + alpha_n*p_n_prev;
        r_n = r_n_prev - alpha_n*A*p_n_prev;
        error(i) = norm(b-A*x_n)^2;
        if norm(r_n) < tol
            return
        end
        beta_n = (transpose(r_n)*r_n)/(transpose(r_n_prev)*r_n_prev);
        p_n = r_n + beta_n*p_n_prev;
    end
end

