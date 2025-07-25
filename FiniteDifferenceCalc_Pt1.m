% Simple egg cooking simulation - Time vs Center Temperature
clear; clc;

% Material properties of Chicken Egg
density = 1175.04; % kg/m³
specific_heat = 3026.05; % J/(kg·K)
k = 0.498; % W/(m·K)
alpha = k/(density * specific_heat); % Thermal diffusivity

% Geometry
r_max = 0.044; % Maximum radius (m)
N = 50; % Spatial grid points (reduced for speed)
dr = r_max / N;

% Time parameters
dt = 0.01; % Time step (s)
t_max = 9000; % Total time (s) - 15 minutes
Nt = round(t_max/dt);

% Spatial grid
r = linspace(dr, r_max, N);

% Initial condition
T0 = 5; % Initial temperature (°C)
T = T0 * ones(1, N);
T_boiling_water = 100; % Surface temperature
T_cooked = 80; % Cooking temperature

% Storage for center temperature vs time
center_temp = zeros(Nt+1, 1);
time_vec = (0:Nt) * dt / 60; % Convert to minutes
center_temp(1) = T(1);



cooked = false;

for n = 1:Nt % Time stepping
    T_old = T;
    
    % Interior points
    for i = 2:N-1
        d2T_dr2 = (T_old(i+1) - 2*T_old(i) + T_old(i-1)) / dr^2;
        dT_dr = (T_old(i+1) - T_old(i-1)) / (2 * dr);
        T(i) = T_old(i) + dt * alpha * (d2T_dr2 + (2/r(i)) * dT_dr);
    end
    
    % Boundary conditions
    T(1) = T_old(1) + dt * alpha * 2 * (T_old(2) - T_old(1)) / dr^2; % Center
    T(N) = T_boiling_water; % Surface
    
    % Store center temperature
    center_temp(n+1) = T(1);
    
    % Check if cooked
    if ~cooked && T(1) >= T_cooked
        cooked = true;
        fprintf('Egg is cooked at %.1f minutes! Center temp: %.1f°C\n', ...
                n*dt/60, T(1));
    end

end

% Plot results
figure;
plot(time_vec, center_temp, 'r-', 'LineWidth', 2);
hold on;
yline(T_cooked, 'k--', 'LineWidth', 1.5, 'Label', 'Cooking Temperature (80°C)');
yline(T0, 'b--', 'LineWidth', 1, 'Label', 'Initial Temperature (5°C)');

xlabel('Time (minutes)');
ylabel('Center Temperature (°C)');
title('Egg Cooking Progress - Center Temperature vs Time');
grid on;
xlim([0, t_max/60]);
ylim([0, 105]);

% Mark cooking point if reached
if cooked
    cook_idx = find(center_temp >= T_cooked, 1);
    cook_time_min = time_vec(cook_idx);
    plot(cook_time_min, T_cooked, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    text(cook_time_min + 0.5, T_cooked + 5, ...
         sprintf('Cooked!\n%.1f min', cook_time_min), ...
         'FontSize', 12, 'Color', 'green');
end

legend('Center Temperature', 'Location', 'southeast');

% Final summary
fprintf('\nSimulation complete!\n');
fprintf('Final center temperature: %.1f°C\n', center_temp(end));
if cooked
    cook_time_final = time_vec(find(center_temp >= T_cooked, 1));
    fprintf('Egg cooked in %.1f minutes\n', cook_time_final);
else
    fprintf('Egg not fully cooked after %.1f minutes\n', t_max/60);
end