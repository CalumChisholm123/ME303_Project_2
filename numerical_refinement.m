%Boundary conditions
BCs = [0,2];
L = 1;
max_time = 10; %seconds
IC = @(x) cos(pi * x);
alpha_sq = 2;

%Fixing delta_t and changing delta_x
delta_t = 0.001;
t_span = 0:delta_t:max_time;
delta_x_values = [0.065, 0.1, 0.2];

figure
for i=1:length(delta_x_values)
    delta_x = delta_x_values(i);
    x_span = 0:delta_x:L;
    %Iterate through rest of the 'grid' of time and space
    u = zeros(length(t_span), length(x_span)); %temp for temperature
    %Apply initial condition
    u(1, :) = IC(x_span);
    %Apply boundary conditions
    u(:,1) = BCs(1);
    u(:,length(x_span)) = BCs(2);
    for k=2:length(t_span)
       for n=2:length(x_span)-1
            u(k, n) = u(k-1, n) + alpha_sq*delta_t/(delta_x^2) * (u(k - 1,n + 1) - alpha_sq*u(k - 1, n) + u(k -1, n - 1));
        end
    end
    plot(x_span, u(t_span == 0.01,:), 'DisplayName', sprintf('\\Deltax = %.3f',delta_x_values(i)))
    hold on
end
xlabel('X position')
ylabel('u(x,t)')
title('Numerical Solution of the 1D Heat Equation - \Deltat = 0.001')
legend
grid on





