%-----------------------------------------------------------------------%
% INTEGRATED MAGNETIC BRAKING SIMULATION
% Combining Coil Field Calculation, Falling Magnet Dynamics, and Induced EMF
%-----------------------------------------------------------------------%
% Rodrigo Gamboa & Francisco Montes | June 2024 (Revised by Gemini AI)
%-----------------------------------------------------------------------%

clear all;          % Clear all variables from the workspace
clc;                % Clear command window
close all;          % Close all figure windows

% --- Coil Geometry and Properties (from Code 1) ---
nl = 1;             % Number of wire loops in the coil
ds_coil = 0.1;      % Length differential for coil discretization (for visualization and B-field calc)
R_coil = 1.5;       % Radius of the conductive ring/coil in meters (used for force and flux calculation)
I_initial_coil = 300; % Initial current for coil visualization (will be replaced by I_induced)
R_resistance = 0.1;   % Resistance of the conductive coil (Ohms) - NEW VARIABLE

% Define 3D working space for coil field visualization
x_space = -5:ds_coil:5;
y_space = -5:ds_coil:5;
z_space = x_space;

Lx_space = length(x_space);
Ly_space = length(y_space);
Lz_space = length(z_space);

rw = 0.2;           % Wire thickness (for visualization, not directly used in force/flux for point dipole)

% Magnetic constants
mo = 4*pi*1e-7;     % Permeability of free space (H/m)

% Parameters for discretizing the coil for visualization
N_coil_points = 100; % Number of points per loop for visualization
sz = 1;             % Loop step size in z-direction (for solenoid stacking)

s_coil_idx = 1;     % Starting index for point arrays

dtheta_coil = 2*pi/N_coil_points; % Angular step for coil discretization
dl_coil = R_coil * dtheta_coil;   % Length of current differential (dl)

ang_coil = 0:dtheta_coil:2*pi-dtheta_coil; % Angles for points around the coil

% --- Calculate and Visualize the Coil (Initial State - No Induced Current Yet) ---
% This part calculates the position vectors for the coil's segments, primarily for visualization.
Px_coil = zeros(1, nl * N_coil_points);
Py_coil = zeros(1, nl * N_coil_points);
Pz_coil = zeros(1, nl * N_coil_points);
dx_coil = zeros(1, nl * N_coil_points);
dy_coil = zeros(1, nl * N_coil_points);
dz_coil = zeros(1, nl * N_coil_points);

for i_loop = 1:nl
    % Calculates x-coordinates of loop points
    Px_coil(s_coil_idx:s_coil_idx+N_coil_points-1)=R_coil*cos(ang_coil);
    % Calculates y-coordinates of loop points
    Py_coil(s_coil_idx:s_coil_idx+N_coil_points-1)=R_coil*sin(ang_coil);
    % Calculates z-coordinates for stacked loops
    Pz_coil(s_coil_idx:s_coil_idx+N_coil_points-1)=-nl/2*sz+(i_loop-1)*sz;

    % Calculates x-component of differential current element (for visualization)
    dx_coil(s_coil_idx:s_coil_idx+N_coil_points-1)=-Py_coil(s_coil_idx:s_coil_idx+N_coil_points-1)*dtheta_coil;
    % Calculates y-component of differential current element (for visualization)
    dy_coil(s_coil_idx:s_coil_idx+N_coil_points-1)=Px_coil(s_coil_idx:s_coil_idx+N_coil_points-1)*dtheta_coil;

    s_coil_idx = s_coil_idx + N_coil_points; % Updates index for next loop
end
dz_coil(1:N_coil_points*nl)=0; % Sets z-component of current differential to zero (planar coil)

% --- Falling Magnet Properties ---
mag = 5e5;     % Magnetic moment of the falling magnet (A*m^2)
m_mass = 0.0365; % Mass of the magnet in kg
g = -9.81;      % Acceleration due to gravity (m/s^2)
w = m_mass * g; % Weight of the magnet in Newtons

zo = 5;         % Initial z-position of the magnet (m)
zring = 0;      % Z-position of the conductive coil (m)

dt = 0.01;      % Time step for the simulation (s)

% Initialize vectors for magnet's motion
t = []; t(1) = 0;       % Time vector
zm = []; zm(1) = zo;    % Magnet's z-position
vz = []; vz(1) = 0;     % Magnet's z-velocity

% Initialize vectors for free-fall comparison
zmfree = []; zmfree(1) = zo;
vzfree = []; vzfree(1) = 0;

% Initialize vectors for induced current and force
I_induced = []; I_induced(1) = 0;
fem = []; fem(1) = 0;
Fm = []; Fm(1) = 0; % Magnetic force on the magnet
F_net = []; F_net(1) = w; % Net force on the magnet (initially just weight)

cc = 1;               % Counter for simulation steps

% --- Setup Figure for Animation ---
figure(1);
clf; % Clear the figure for animation
hold on;
grid on;
xlabel('X Position (m)');
ylabel('Z Position (m)');
title('Magnetic Braking Simulation');
ylim([-R_coil*2 zo + 1]); % Adjust Y-axis limits dynamically based on R_coil and zo
xlim([-R_coil*2 R_coil*2]); % Adjust X-axis limits based on R_coil

% Plot the stationary coil
plot(Px_coil, Pz_coil, 'k', 'LineWidth', 2); % Plotting the coil structure
quiver3(Px_coil(1:10:end), Py_coil(1:10:end), Pz_coil(1:10:end), ...
        dx_coil(1:10:end), dy_coil(1:10:end), dz_coil(1:10:end), ...
        0.5, 'k'); % Arrows indicating coil's orientation (optional, if I_induced direction is fixed)

% Animated lines for paths
path_braked = animatedline('Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Braked Fall');
path_free = animatedline('Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Free Fall');

% Legend for animated lines
legend('show', 'Location', 'northwest');

% --- Simulation Loop ---
% The magnet falls until it reaches just before the coil center (z=0)
while (zm(cc) > -10) % Stop point as specified (just before passing through coil)

    % 1. Calculate magnetic flux at current and next time step (due to falling magnet)
    % Call B_due_M for current position
    [~,~,phiB1,~] = B_due_M(zm(cc), mag, R_coil);

    % 2. Update magnet's position and velocity (temporarily for free fall)
    % These are initial kinematic updates assuming only gravity for phiB2 calculation
    % The force from induced current will be applied after fem calculation
    temp_vz_next = vz(cc) + g*dt; % Simple Euler for velocity update (vz + a*dt)
    temp_zm_next = zm(cc) + vz(cc)*dt + 0.5*g*dt^2; % Simple Euler for position update (z + v*dt + 0.5*a*dt^2)

    % Ensure temp_zm_next doesn't go below the stop threshold if dt is too large
    %if temp_zm_next < 0.0162
    %    temp_zm_next = 0.0162;
    %end

    % Call B_due_M for next estimated position
    [~,~,phiB2,~] = B_due_M(temp_zm_next, mag, R_coil);

    % 3. Calculate Induced EMF (Faraday's Law)
    fem(cc) = (phiB2 - phiB1)/dt; % Negative sign for Lenz's Law 
                                  % Removed negative due nonsense.

    % 4. Calculate Induced Current
    I_induced(cc) = fem(cc) / R_resistance;

    % 5. Calculate Magnetic Force on the Falling Magnet (due to I_induced in coil)
    % Force between a magnetic dipole and a current loop
    % Note: I_induced can be negative, which is correct for direction.
    Fm(cc) = (3 * zm(cc) * mag * mo * I_induced(cc) * R_coil^2) / (2 * (zm(cc)^2 + R_coil^2)^(5/2));

    % 6. Calculate Net Force and Acceleration
    F_net(cc) = Fm(cc) + w; % Net force = Magnetic Force + Weight
    a_braked = F_net(cc) / m_mass;

    % 7. Update Magnet's Position and Velocity (with magnetic braking)
    vz(cc+1) = vz(cc) + a_braked * dt;
    zm(cc+1) = zm(cc) + vz(cc) * dt + 0.5 * a_braked * dt^2;

    % 8. Update Free Fall Position and Velocity (for comparison)
    vzfree(cc+1) = vzfree(cc) + g * dt;
    zmfree(cc+1) = zmfree(cc) + vzfree(cc) * dt + 0.5 * g * dt^2;

    % Ensure free fall doesn't go past the target if the main sim stops
    %if zmfree(cc+1) < 0.0162
    %     zmfree(cc+1) = 0.0162;
    %     vzfree(cc+1) = vzfree(cc); % Stop its velocity too
    %end

    % 9. Update Time
    t(cc+1) = t(cc) + dt;

    % 10. Animation and Display
    fprintf('Time: %.2f s, Z position: %.4f m, Induced Current: %.4e A\n', t(cc), zm(cc), I_induced(cc));

    % Plot current magnet positions
    h_braked = scatter(0, zm(cc), 150, 'b', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Braked Magnet');
    h_free = scatter(0, zmfree(cc), 150, 'r', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', 'Free Fall Magnet');

    % Add points to animated lines
    addpoints(path_braked, 0, zm(cc));
    addpoints(path_free, 0, zmfree(cc));

    drawnow limitrate; % Update plot efficiently

    % Delete current scatter points for next frame
    delete(h_braked);
    delete(h_free);

    cc = cc + 1; % Increment step counter
end

hold off; % Release the plot hold

% --- Post-Simulation Plotting ---
figure(2); % Create new figure for results
sgtitle('Magnetic Braking Simulation Results'); % Super title for the figure

% Plot 1: Magnetic Force vs. Z Position
subplot(2,3,1);
plot(zm(1:length(Fm)), 1000*Fm, '-b', 'LineWidth', 2);
hold on;
plot([max(zm), min(zm)], [0,0], '-.k', 'LineWidth', 1); % Line at Fm=0
grid on;
xlabel('Z Position (m)');
ylabel('Magnetic Force (mN)');
title('Magnetic Force on Falling Magnet');

% Plot 2: Net Force vs. Z Position
subplot(2,3,2);
plot(zm(1:length(F_net)), 1000*F_net, '-b', 'LineWidth', 2);
hold on;
plot([max(zm), min(zm)], [0,0], '-.k', 'LineWidth', 1); % Line at Fm=0
grid on;
xlabel('Z Position (m)');
ylabel('Net Force (mN)');
title('Net Force on Falling Magnet');

% Plot 3: Induced EMF vs. Time
subplot(2,3,3);
plot(t(1:length(fem)), fem, '-g', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Induced EMF (V)');
title('Induced Electromotive Force');

% Plot 4: Z Position vs. Time
subplot(2,3,4);
hold on;
plot(t, zm, '-b', 'LineWidth', 2, 'DisplayName', 'Braked Fall');
plot(t, zmfree, '--r', 'LineWidth', 2, 'DisplayName', 'Free Fall');
plot([0, max(t)], [zring, zring], '-.k', 'LineWidth', 1, 'DisplayName', 'Coil Location');
grid on;
xlabel('Time (s)');
ylabel('Z Position (m)');
title('Position vs. Time');
legend('show', 'Location', 'southwest');

% Plot 5: Induced Current vs. Time
subplot(2,3,5);
plot(t(1:length(I_induced)), I_induced, '-m', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Induced Current (A)');
title('Induced Current in Coil');


%--------------------------------------------------------------------------
% Function to calculate magnetic field (Bz) and flux (phiB) due to a
% point magnetic dipole (magnet) through a circular area (coil)
%--------------------------------------------------------------------------
function [x_grid, y_grid, phiB, Bz_field]=B_due_M(z_magnet, mag_dipole, R_ring)

    ds_grid = 0.005; % Grid step size for numerical integration of flux
    
    x_grid = -R_ring:ds_grid:R_ring; % X-coordinates for the grid over the coil's area
    y_grid = x_grid;                % Y-coordinates for the grid over the coil's area

    Lx_grid = length(x_grid);
    Ly_grid = length(y_grid);

    Bz_field = zeros(Lx_grid, Ly_grid); % Initialize 2D matrix for z-component of B-field
    phiB = 0; % Initialize magnetic flux

    mo = 4*pi*1e-7; % Permeability of free space

    for i_x = 1:Lx_grid
        for j_y = 1:Ly_grid
            % Radial distance from the center of the coil in the x-y plane
            r_xy = sqrt(x_grid(i_x)^2 + y_grid(j_y)^2);

            % Only consider points within the ring's area
            if r_xy <= R_ring

                % Distance from the dipole (at (0,0,z_magnet)) to the current point on the coil plane (x_grid(i_x), y_grid(j_y), 0)
                % r_vector = (x_grid(i_x) - 0)i + (y_grid(j_y) - 0)j + (0 - z_magnet)k
                % Components of the position vector from dipole to field point
                rx_from_dipole = x_grid(i_x);
                ry_from_dipole = y_grid(j_y);
                rz_from_dipole = -z_magnet; % The z-coordinate of the field point is 0 (coil plane) relative to magnet at z_magnet

                % Magnitude of the distance vector from dipole to field point
                r_mag_from_dipole = sqrt(rx_from_dipole^2 + ry_from_dipole^2 + rz_from_dipole^2);

                % The original formula from function for Bz is:
                % (mo/(4*pi)) * ((3*z*(mag*z)-mag*(x(i)^2 + y(j)^2 + z^2))/(x(i)^2 + y(j)^2 + z^2 + (ds^2))^(-5/2));
                % This translates to (mo/(4*pi)) * (3*z_magnet^2 - (r_xy^2 + z_magnet^2)) / (r_xy^2 + z_magnet^2 + ds_grid^2)^(5/2)
                Bz_field(i_x, j_y) = (mo/(4*pi)) * ...
                                     ( (3*z_magnet*(mag_dipole*z_magnet) - mag_dipole*(x_grid(i_x)^2 + y_grid(j_y)^2 + z_magnet^2)) / ...
                                       (x_grid(i_x)^2 + y_grid(j_y)^2 + z_magnet^2 + (ds_grid^2))^(5/2) );

                % Accumulate magnetic flux using numerical integration (Area * B_z)
                phiB = phiB + ds_grid^2 * Bz_field(i_x, j_y);
            end
        end
    end
end