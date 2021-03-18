close all
%generate matrix A
arr_A = cell(3,2);
m_vals = [100, 500, 2500];
tau_vals = [0.1, 0.01];
for m_i = 1 : size(m_vals, 2)
    for tau_j = 1 : size(tau_vals, 2)
        arr_A{m_i, tau_j} = produce_S(m_vals(m_i), tau_vals(tau_j));
    end
end
sigma_w_values = [0.0001, 0.01, 1];
arr_b = cell(size(sigma_w_values, 2),10);
arr_x_0 = cell(10,1);
%%
E_S_psuedo = zeros(6, 3);
E_S_CG = zeros(6, 3);
E_S_GMRES = zeros(6, 3);
E_S_MINRES = zeros(6, 3);
for A_i = 1:6
    E_S_sum_psuedo = zeros(1, 3);
    E_S_sum_CG = zeros(1, 3);
    E_S_sum_GMRES = zeros(1, 3);
    E_S_sum_MINRES = zeros(1, 3);
    A = arr_A{idivide(int32(A_i), int32(2),'ceil'),2-rem(A_i,2)};
    %
    m = m_vals(idivide(int32(A_i), int32(2),'ceil'));
    for i = 1:10
        x_0 = randn(m, 1);
        b_0 = A*x_0;
        for sigma_w_j = 1:size(sigma_w_values, 2)
            sigma_w = sigma_w_values(sigma_w_j);
            w = sigma_w * randn(m, 1);
            arr_b{sigma_w_j, i} = b_0 + w;
        end
        arr_x_0{i} = x_0;
    end
    %
    for i = 1:3
        for j = 1:10
            x_0_j = arr_x_0{j};
            b_i_j = arr_b{i, j};
            x_hat_i_j_psuedo = psuedo_inverse(A)*b_i_j;
            x_hat_i_j_CG = CG(A, b_i_j, 25, eps);
            x_hat_i_j_GMRES = GMRES(A, b_i_j, 25, eps);
            x_hat_i_j_MINRES = MINRES(A, b_i_j, 25, eps);
            %error summation
            E_S_sum_psuedo(i) = E_S_sum_psuedo(i) + norm(-x_hat_i_j_psuedo + x_0_j)^2;
            E_S_sum_CG(i) = E_S_sum_CG(i) + norm(-x_hat_i_j_CG + x_0_j)^2;
            E_S_sum_GMRES(i) = E_S_sum_GMRES(i) + norm(-x_hat_i_j_GMRES + x_0_j)^2;
            E_S_sum_MINRES(i) = E_S_sum_MINRES(i) + norm(-x_hat_i_j_MINRES + x_0_j)^2;
            i
            j
            A_i
        end
    end
    E_S_psuedo(A_i,:) = sqrt( 0.1 * E_S_sum_psuedo);
    E_S_CG(A_i,:) = sqrt( 0.1 * E_S_sum_CG);
    E_S_GMRES(A_i,:) = sqrt( 0.1 * E_S_sum_GMRES);
    E_S_MINRES(A_i,:) = sqrt( 0.1 * E_S_sum_MINRES);
end
%plot E_S_psuedo
figure
title('E_{S,i} with respect to \sigma_w for Pseudoinverse');
ylabel('E_{S,i}');
xlabel('\sigma_w');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold on
for i = 1:6
    plot(sigma_w_values, E_S_psuedo(i,:));
end
legend('m=100, \tau=0.1','m=100, \tau=0.01', 'm=500, \tau=0.1', 'm=500, \tau=0.01', 'm=2500, \tau=0.1', 'm=2500, \tau=0.01')
hold off
%plot E_S_CG
figure
title('E_{S,i} with respect to \sigma_w for CG');
ylabel('E_{S,i}');
xlabel('\sigma_w');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold on
for i = 1:6
    plot(sigma_w_values, E_S_CG(i,:));
end
legend('m=100, \tau=0.1','m=100, \tau=0.01', 'm=500, \tau=0.1', 'm=500, \tau=0.01', 'm=2500, \tau=0.1', 'm=2500, \tau=0.01')
hold off
%plot E_S_GMRES
figure
title('E_{S,i} with respect to \sigma_w for GMRES');
ylabel('E_{S,i}');
xlabel('\sigma_w');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold on
for i = 1:6
    plot(sigma_w_values, E_S_GMRES(i,:));
end
legend('m=100, \tau=0.1','m=100, \tau=0.01', 'm=500, \tau=0.1', 'm=500, \tau=0.01', 'm=2500, \tau=0.1', 'm=2500, \tau=0.01')
hold off
%plot E_S_MINRES
figure
title('E_{S,i} with respect to \sigma_w for MINRES');
ylabel('E_{S,i}');
xlabel('\sigma_w');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold on
for i = 1:6
    plot(sigma_w_values, E_S_MINRES(i,:));
end
legend('m=100, \tau=0.1','m=100, \tau=0.01', 'm=500, \tau=0.1', 'm=500, \tau=0.01', 'm=2500, \tau=0.1', 'm=2500, \tau=0.01')
hold off
%%
arr_E_O = cell(6,1);
max_iter_no = 20;
for A_i = 1:6
    A = arr_A{idivide(int32(A_i), int32(2),'ceil'),2-rem(A_i,2)};
    %
    m = m_vals(idivide(int32(A_i), int32(2),'ceil'));
    for i = 1:10
        x_0 = randn(m, 1);
        b_0 = A*x_0;
        for sigma_w_j = 1:size(sigma_w_values, 2)
            sigma_w = sigma_w_values(sigma_w_j);
            w = sigma_w * randn(m, 1);
            arr_b{sigma_w_j, i} = b_0 + w;
        end
        arr_x_0{i} = x_0;
    end
    %
    E_O = zeros(3, max_iter_no);
    for i = 1:3
        for j = 1:10
            b_i_j = arr_b{i, j};
            [x, iter, error] = CG(A, b_i_j, max_iter_no, eps);
            E_O(i,:) = E_O(i,:) + error;
        end
        E_O(i,:) = sqrt(0.1 * E_O(i,:));
    end
    arr_E_O{A_i} = E_O;
end
for m_i=1:3
    for sigma_i=1:3
        figure
        sigma = sigma_w_values(sigma_i);
        m = m_vals(m_i);
        title(['E_{O,i} with respect to n for CG where m=',num2str(m),' and \sigma=', num2str(sigma)]);
        ylabel('E_{O,i}');
        xlabel('n');
        set(gca, 'YScale', 'log')
        hold on
        for tau_i=0:1
            E_O = arr_E_O{2*m_i-1+tau_i};
            plot((1:max_iter_no), E_O(sigma_i,:));
        end
        legend('\tau=0.1','\tau=0.01')
        hold off
    end
end
%%
% b = arr_b{2,1};
% %Pseudoinverse
% A_psuedo = psuedo_inverse(A);
% x = A_psuedo*b;
% %CG
% [x, i, error] = CG(A, b, 20, eps);
% figure
% plot((1:i), error(1:i));
% title('Convergence Curve for CG');
% ylabel('||r_n||');
% xlabel('n');
% set(gca, 'YScale', 'log')
% %GMRES
% [x, i, error] = GMRES(A, b, 20, eps);
% figure
% plot([1:i], error(1:i));
% title('Convergence Curve for GMRES');
% ylabel('||r_n||');
% xlabel('n');
% set(gca, 'YScale', 'log')
% %MINRES
% [x, i, error] = MINRES(A, b, 20, eps);
% figure
% plot([1:i], error(1:i));
% title('Convergence Curve for MINRES');
% ylabel('||r_n||');
% xlabel('n');
% set(gca, 'YScale', 'log')
