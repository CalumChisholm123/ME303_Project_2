% Simple egg cooking simulation - Time vs Center Temperature
clear; clc; close;

% Material properties of Egg
density = 1150; % kg/m³
specific_heat = 3397; % J/(kg·K)
k = 0.535; % W/(m·K)
alpha = k/(density * specific_heat); % Thermal diffusivity

% Geometry
r_max = 0.12954/2; % Maximum radius (m)
N = 50; % Spatial grid points (reduced for speed)
dr = r_max / N;

% Time parameters
dt = 0.01; % Time step (s)
t_max = 9000; % Total time 
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
time_vec = (0:Nt) * dt / 60; % 
center_temp(1) = T(1);

cooked = false;
cooked_duration = 0;  % Time egg is cooked 
required_cooked_time = 10;  % 10 seconds

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
    
    center_temp(n+1) = T(1);
    
    % Check if cooked 
    if T(1) >= T_cooked
        if ~cooked
            % Reached 80 degrees Celcius 
            cooked = true;
            cooked_start_time = n * dt;
        end

        cooked_duration = (n * dt) - cooked_start_time;
    else
        % Temperature dropped below threshold, reset cooking
        cooked = false;
        cooked_duration = 0;
    end
    
    % Check if cooked for required time
    if cooked && cooked_duration >= required_cooked_time
        fprintf('Egg has been cooked. Center temp: %.1f°C\n', T(1));
    end
end

if cooked_duration < required_cooked_time
    fprintf('Egg was not cooked for at least 10 seconds');
end

% Plotting
figure;
plot(time_vec, center_temp, 'r-', 'LineWidth', 1);
hold on;
yline(T_cooked, 'k--', 'LineWidth', 1);
yline(T0, 'b--', 'LineWidth', 1);

xlabel('Time (minutes)');
ylabel('Center Temperature (°C)');
title('Ostrich Egg Center Temperature vs Time');
grid on;
xlim([0, t_max/60]);
ylim([0, 105]);

legend('Center Temperature', 'Location', 'northwest');

% Cooked output 
fprintf('Final  temperature of center of egg is: %.1f°C\n', center_temp(end));
if cooked
    cook_time_final = time_vec(find(center_temp >= T_cooked, 1));
    fprintf('Egg cooked in %.1f minutes\n', cook_time_final);
else
    fprintf('Egg was not cooked after %.1f minutes\n', t_max/60);
end