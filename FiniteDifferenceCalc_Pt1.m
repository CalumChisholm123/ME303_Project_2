
density = 1009.62; %kg/m^3
specific_heat = 498; %J/kgK
k = 3.163; % W/(mK)
alpha = k/(density * specific_heat);       % Thermal diffusivity

fprintf('Thermal diffusivity: %.2e m²/s\n', alpha);

r_max = 0.044;       % Maximum radius
N = 100;           % Spatial grid points
dr = r_max / N;    % Spatial step
dt = 0.0001;       % Time step
t_max = 600;       % Total time

% Spatial grid (avoid r=0)
r = linspace(dr, r_max, N);  

% Initial condition
T0 = 5;  
T = T0 * ones(1, N); 

cooked = false;
cook_time = 0;
T_cooked = 80;

M = length(r);
T_new = zeros(1, M);

for t = 0:dt:t_max % Time loop 
    
    for i = 2:M-1
        d2T_dr2 = (T(i+1) - 2*T(i) + T(i-1)) / dr^2; % 2nd order central difference 
        dT_dr = (T(i+1) - T(i-1)) / (2 * dr); % 1st order central difference 
        T_new(i) = T(i) + dt * alpha * (d2T_dr2 + (2/r(i)) * dT_dr); % Computing new Temperature 
        
    % Check if egg is cooked (center reaches 80°C)
    if ~cooked && T_new(1) >= T_cooked
        cooked = true;
        cook_time = t * dt;
        fprintf('Egg is cooked! Time: %.1f seconds (%.1f minutes)\n', ...
                cook_time, cook_time/60);
        fprintf('Center temperature: %.1f°C\n', T_new(1));
    end


    end
    
    % Boundary conditions
    T_new(1) = T_new(2);       % Symmetry at r=0 (Neumann)
    T_new(M) = 100;              % Fixed temperature at r=r_max (Dirichlet)
    
    T = T_new;
end

% Plot
Initial_Profile = T0 * ones(1, N);
plot(r, Initial_Profile, 'k--', r, T, 'r-', 'LineWidth', 2);
xlabel('Radius (r)');
ylabel('Temperature (T)');
legend('Initial', 'Final');
title('Forward Euler Solution of Spherical Heat Equation');
grid on;