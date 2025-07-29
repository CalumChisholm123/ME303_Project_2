
clc;            
clear;          
close all;      

% Define the spatial domain for x
x = linspace(0, 1, 500); % Vector of 500 points from x=0 to x=1

% Define the time values from the legend
t_values = [0.001, 0.01, 0.1, 10];

N = 150; 

% --- Plotting ---
figure; % Create a new figure window
hold on; % Hold the plot to overlay multiple lines

% Loop through each specified time value
for i = 1:length(t_values)
    t = t_values(i);
    
    % Initialize the transient part of the solution, v(x,t)
    v = zeros(size(x)); 
    
    % Calculate the series summation for v(x,t)
    for n = 1:N
        % Calculate the coefficient C_n based on its piecewise definition
        if n == 1
            Cn = -4/pi;
        else % n > 1
            % NOTE: Corrected formula is used here. See explanation below.
            term1 = (2 * n * ((-1)^n + 1)) / (pi * (n^2 - 1));
            term2 = -4 * ((-1)^(n+1)) / (n * pi);
            Cn = term1 + term2;
        end
        
        % Add the current term to the transient solution v
        v = v + Cn * sin(n*pi*x) .* exp(-2*(n*pi)^2 * t);
    end
    
    % Calculate the full solution u(x,t) = v(x,t) + u_E(x)
    uE = 2*x; % Steady-state part of the solution
    u = v + uE;
    
    % Plot the result for the current time t
    plot(x, u, 'LineWidth', 2, 'DisplayName', sprintf('t = %gs', t));
end

% --- Finalize the Plot ---
hold off; % Release the plot hold
grid on; % Add a grid
xlabel('x'); % Label the x-axis
ylabel('u(x,t)'); % Label the y-axis
title('Plot of u(x,t) vs. x for Different Timesteps'); % Add a title
legend('show', 'Location', 'northwest');