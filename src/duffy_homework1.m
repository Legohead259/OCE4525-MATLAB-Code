% OCE4525 - Homework 2 - Wave Transformations
% Braidan Duffy
% 09/08/2021

clear all; close all;

% Perform the calculations and place the answers in the table below so you can compare what the ocean wave parameters 
% would be in 3 different depths of 200m, 10m, and 2m.  Use linear wave theory for a 1.5 meter high wave (H = 1.5m) with a
% 12-second period (T = 12sec).
% Use:  ρ = 1025 kg/m3 /  g = 9.81 m/s2

d = [200, 10, 2]; % m
H = 1.5; % m
T = 12; % s
rho_sw = 1025; % kg/m^3
g = 9.81; % m/s^2
len = length(d);

%% Part A - General wave parameters
col_names = {'d=200m', 'd=10m', 'd=2m'};
row_names = {'Wavelength [m]', 'Relative Depth', 'Celerity [m/s]', 'Group Celerity [m/s]', 'Energy [J]', 'Power [W]', 'Hmax [m]', 'Ursell Number'};
parta_table = table('RowNames', row_names);

for i=1:len % For every depth
    [L, C, ~, Cg] = wave_params(H, T, d(i), 0, rho_sw, g);    % Calculate wavelength. celerity, group celerity, wave energy, wave power
    rel_depth = d(i)/L;                                       % Calculate relative depth
    E = 1/8 * rho_sw * g * H^2 * L * 100;                     % Calculate wave energy for 1 wave with a 100m crest width
    F = E * Cg / L;                                           % Calculate wave power for 1 wave with a 100m crest width
    Hmax_depth = 0.78 * d(i);                                 % Calculate maximum wave height due to depth
    Hmax_steep = L / 7;                                       % Calculate maximum wave height due to wave steepness
    Hmax = min(Hmax_depth, Hmax_steep);                       % Determine maximum wave height possible
    Ur = H * L^2 / d(i)^3;                                    % Calculate Ursell number

    parta_table = addvars(parta_table, [L; rel_depth; C; Cg; E; F; Hmax; Ur], 'NewVariableNames', col_names(i));
end

disp(parta_table)

%% Part B - Particle and Pressures at 2m below SWL
z = -2.0; % m
row_names = {'Hydrostatic Pressure [Pa]', 'Max Dynamic Pressure [Pa]', 'Max Gauge Pressure [Pa] (crest)', 'Max Gauge Pressure [Pa] (trough)'...
                'Max Horizontal Particle Velocity [m/s]', 'Max Vertical Particle Velocity [m/s]', ...
                'Max Horizontal Particle Acceleration [m/s/s]', 'Max Vertical Particle Acceleration [m/s/s]', ...
                'Max Horizontal Particle Displacement [m]', 'Max Vertical Particle Displacement [m]'};
partb_table = table('RowNames', row_names);

for i=1:len % For every depth
    [Sx, Sz, u, w, ax, az, ~, Pd, Ph, P_crest, ~] = wave_particles(H, T, d(i), z, g, rho_sw);   % Calculate particle motions and pressures
    P_trough = Ph - Pd;                                                                         % Calculate maximum gauge pressure beneath trough
    partb_table = addvars(partb_table, [Ph; Pd; P_crest; P_trough; u; w; ax; az; Sx; Sz], 'NewVariableNames', col_names(i)); % Add data to table
end

disp(partb_table)