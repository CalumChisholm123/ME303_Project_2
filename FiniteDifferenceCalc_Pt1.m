% Corrected egg cooking simulation using spherical heat equation
clear; clc;

% Material properties (corrected units)
density = 1009.62; % kg/m³
specific_heat = 498; % J/(kg·K) - converted from kJ to J
k = 3.163; % W/(m·K)
alpha = k/(density * specific_heat); % Thermal diffusivity [m²/s]

fprintf('Thermal diffusivity: %.2e m²/s\n', alpha);

% Geometry
r_max = 0.044; % Maximum radius (m) - about 4.4 cm
N = 100; % Spatial grid points
dr = r_max / N; % Spatial step

% Time parameters - much longer simulation time
dt = 0.01; % Time step (s) - adjusted for stability
t_max = 600; % Total time (s) - 10 minutes
Nt = round(t_max/dt);

% Check stability (Courant number should be < 0.5 for stability)
courant = alpha * dt / dr^2;
fprintf('Courant number: %.4f (should be < 0.5)\n', courant);

if courant >= 0.5
    warning('Time step may be too large for stability!');
    dt = 0.4 * dr^2 / alpha; % Adjust for stability
    fprintf('Adjusted time step: %.4f s\n', dt);
    Nt = round(t_max/dt);
end

% Spatial grid (avoid r=0)
r = linspace(dr, r_max, N);

% Initial condition
T0 = 5; % Initial temperature (°C)
T = T0 * ones(1, N);
T_boiling_water = 100; % Boiling water temperature (°C)
T_cooked = 80; % Temperature when egg is considered cooked (°C)

% Storage for animation/plotting
T_history = zeros(Nt+1, N);
T_history(1, :) = T;
time_vec = 0:dt:(Nt*dt);

fprintf('Starting simulation...\n');

% Time stepping
cooked = false;
cook_time = 0;

fprintf('Initial center temperature: %.1f°C\n', T(1));

for n = 1:Nt
    T_old = T; % Store current temperature
    
    % Update interior points (i = 2 to N-1)
    for i = 2:N-1
        d2T_dr2 = (T_old(i+1) - 2*T_old(i) + T_old(i-1)) / dr^2;
        dT_dr = (T_old(i+1) - T_old(i-1)) / (2 * dr);
        T(i) = T_old(i) + dt * alpha * (d2T_dr2 + (2/r(i)) * dT_dr);
    end
    
    % Boundary condition at center (i=1): symmetry dT/dr = 0
    % Use the heat equation at r = dr with the symmetry condition
    % At r=dr, dT/dr ≈ (T(2) - T(0))/2dr = (T(2) - T(1))/2dr = 0 (by symmetry)
    % This gives us T(0) = T(2), but we're at i=1 which represents r=dr
    % So we use: d²T/dr² = (T(2) - 2*T(1) + T(0))/dr² = (T(2) - 2*T(1) + T(2))/dr² = 2*(T(2) - T(1))/dr²
    d2T_dr2_center = 2 * (T_old(2) - T_old(1)) / dr^2;
    T(1) = T_old(1) + dt * alpha * d2T_dr2_center;
    
    % Boundary condition at surface (i=N): fixed temperature
    T(N) = T_boiling_water;
    
    % Check if egg is cooked (center reaches 80°C)
    if ~cooked && T(1) >= T_cooked
        cooked = true;
        cook_time = n * dt;
        fprintf('Egg is cooked! Time: %.1f seconds (%.1f minutes)\n', ...
                cook_time, cook_time/60);
        fprintf('Center temperature: %.1f°C\n', T(1));
    end
    
    % Store temperature history
    T_history(n+1, :) = T;
    
    % Progress update
    if mod(n, round(Nt/10)) == 0
        fprintf('Progress: %.0f%%, Center temp: %.1f°C\n', ...
                100*n/Nt, T(1));
    end
end

% Final results
fprintf('\nSimulation complete!\n');
fprintf('Final center temperature: %.1f°C\n', T(1));
if cooked
    fprintf('Egg cooked in %.1f minutes\n', cook_time/60);
else
    fprintf('Egg not fully cooked after %.1f minutes\n', t_max/60);
end

% Plotting
figure(1);
subplot(2,1,1);
plot(r*1000, T_history(1,:), 'b--', 'LineWidth', 2); hold on;
plot(r*1000, T_history(end,:), 'r-', 'LineWidth', 2);
if cooked
    cook_idx = round(cook_time/dt) + 1;
    plot(r*1000, T_history(cook_idx,:), 'g-', 'LineWidth', 2);
    legend('Initial', 'Final', sprintf('Cooked (%.1f min)', cook_time/60), 'Location', 'best');
else
    legend('Initial', 'Final', 'Location', 'best');
end
xlabel('Radius (mm)');
ylabel('Temperature (°C)');
title('Temperature Profile in Spherical Egg');
grid on;
ylim([0 105]);

subplot(2,1,2);
plot(time_vec/60, T_history(:,1), 'r-', 'LineWidth', 2); hold on;
plot(time_vec/60, T_history(:,round(N/2)), 'g-', 'LineWidth', 2);
plot(time_vec/60, T_history(:,end), 'b-', 'LineWidth', 2);
if cooked
    plot(cook_time/60, T_cooked, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
xlabel('Time (minutes)');
ylabel('Temperature (°C)');
legend('Center', 'Midpoint', 'Surface', 'Cooked point', 'Location', 'best');
title('Temperature vs Time at Different Positions');
grid on;
yline(T_cooked, 'k--', 'Cooking Temperature');

% Animation (optional)
create_animation = false; % Set to true if you want animation
if create_animation
    figure(2);
    for i = 1:10:length(time_vec)
        plot(r*1000, T_history(i,:), 'r-', 'LineWidth', 2);
        xlabel('Radius (mm)');
        ylabel('Temperature (°C)');
        title(sprintf('Temperature Profile at t = %.1f min', time_vec(i)/60));
        ylim([0 105]);
        grid on;
        pause(0.1);
    end
end