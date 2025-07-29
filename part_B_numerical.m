%Boundary conditions
BCs = [0,2];
L = 1;
max_time = 10; %seconds
IC = @(x) cos(pi * x);
alpha_sq = 2;
delta_t = 0.00001;
delta_x = 0.01;

[u, x_span, t_span] = numerical_solver(BCs, L, max_time, IC, alpha_sq, delta_t, delta_x);

%Find indices of time = 0.001, 0.01, 0.1, 10
figure
plot(x_span, u(t_span == 0.001,:), 'DisplayName', 't = 0.001s')
hold on
plot(x_span, u(t_span == 0.01,:), 'DisplayName', 't = 0.01s')
hold on
plot(x_span, u(t_span == 0.1,:), 'DisplayName', 't = 0.1s')
hold on
plot(x_span, u(t_span == 10,:), 'DisplayName', 't = 10s')
hold on
xlabel('X position')
ylabel('u(x,t)')
title('Numerical Solution of the 1D Heat Equation')
legend('Location', 'northwest')
grid on

%Function - outputs temperature matrix t rows, x columns
function [u, x_span, t_span] = numerical_solver(BCs, bar_length, max_t, IC_function, alpha, delta_t, delta_x)
    x_span = 0:delta_x:bar_length;
    t_span = 0:delta_t:max_t;

    u = zeros(length(t_span), length(x_span)); %temp for temperature
    %Apply initial condition
    u(1, :) = IC_function(x_span);
    %Apply boundary conditions
    u(:,1) = BCs(1);
    u(:,length(x_span)) = BCs(2);
    
    %Iterate through rest of the 'grid' of time and space
    for k=2:length(t_span)
        for n=2:length(x_span)-1
            u(k, n) = u(k-1, n) + alpha*delta_t/(delta_x^2) * (u(k - 1,n + 1) - alpha*u(k - 1, n) + u(k -1, n - 1));
        end
    end
    
end