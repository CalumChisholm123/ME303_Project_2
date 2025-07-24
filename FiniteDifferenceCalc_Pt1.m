%% Finite Difference Method to solve the egg
clc;clear;close

%% Parameters 
R = 0.0225; % Radius of Egg
T_final = 200; %Total time in seconds 
k = 3.163;
density = 1070.9;
specific_heat = 0.498;

alpha = k/(density*specific_heat); % Thermal diffusivity constant 
T_inital =5; % Inital temperature of the egg 
T_surface = 100; % Surface temperature of the egg in boiling water 

%% Discretization 
Nr = 20; % Number of Spatial grid points 
Nt = 5000; % Number of time steps
dr = R/Nr; %Spatial step size 
dt = T_final/Nt; %Time step size 
r = linspace(0,R,Nr+1); % Radial grid 
t = linspace(0,T_final,Nt+1); % Time grid 

%% Stability Check
stability_factor = alpha * dt / dr^2;
fprintf('Stability factor is', stability_factor);
if stability_factor > 0.5
    warning('Solution may be unstable.');
end

%% Initalize Solution Matrix 
U = zeros(Nr+1,Nt+1);
U(1:Nr,1) = T_inital; % Inital conditions
U(Nr+1,:) = T_surface; % Boundary condition 

for n = 1:Nt % Loop through time
    for i = 2:Nr 
        d2U_dr2 = (U(i+1, n) - 2*U(i, n) + U(i-1, n)) / dr^2;
        dU_dr = (2 / r(i)) * (U(i+1, n) - U(i-1, n)) / (2 * dr);
        
        U(i, n+1) = U(i, n) + alpha * dt * (d2U_dr2 + dU_dr);
    end
end


%% Check for "Cooked Egg" Criterion
T_center_history = U(1, :); 
T_cooked = 80; % Your criterion for a cooked yolk

cooked_time_index = find(T_center_history >= T_cooked, 1, 'first');

if ~isempty(cooked_time_index)
    cooked_time_seconds = t(cooked_time_index);
    fprintf('Success! The egg center reaches %d°C at %.1f seconds (%.2f minutes).\n', ...
            T_cooked, cooked_time_seconds, cooked_time_seconds/60);
else
    fprintf('The egg center did not reach %d°C within the simulated time of %.1f minutes.\n', ...
            T_cooked, T_final/60);
end

%% Plotting the Results

% 3D Surface Plot
figure;
surf(r, t, U', 'EdgeColor', 'none');
colorbar;
xlabel('Radius (r) [m]');
ylabel('Time (t) [s]');
zlabel('Temperature (°C)');
title('Temperature Distribution in Egg vs. Time');
view(60, 30); % Adjust view angle

% 2D Plot of Temperature Profiles
figure;
hold on;
% Select 6 evenly spaced time points to plot
plot_indices = round(linspace(1, Nt+1, 6)); 
for idx = plot_indices
    plot(r, U(:, idx), 'LineWidth', 1.5, 'DisplayName', sprintf('t = %.0f s', t(idx)));
end
hold off;
grid on;
xlabel('Radius (r) [m]');
ylabel('Temperature (°C)');
title('Temperature Profile at Different Times');
legend;