% OCE4525 - Homework 2 - Wave Transformations
% Braidan Duffy
% 09/06/2021

clear all; close all;

%% Part A
% Wave Refraction, Shoaling & Breaking: For the design of a structure in water depth = 6m, 
% you need to analyze the possible design waves. Consider a 3-meter high wave in deep-water (Ho = 3.0m ) 
% with periods of 8 and 16 seconds (T1=8 sec and T2 = 16 sec) with three different angles to the straight 
% parallel offshore depth contours of 0, 30 and 45 degrees, where the offshore slope is 1 vertical to 20 horizontal 
% (1V:20H or m = tan β = 1/20 = 0.05 or cot β = 20):

d = 6.0;                % water depth [m]
H0 = 3.0;               % deepwater wave height [m]
periods = [8, 16];      % wave periods [s]
alphas = [0, 30, 45];   % deep water wave angles [deg]
m = 1/20;               % beach slope
beta = atand(m);        % beach slope [deg]

super_col_names = {'α=0°', 'α=30°', 'α=45°'};
col_names = {'T=8s', 'T=16s'};
row_names = {'Wavelength [m]', 'Unrefracted H0 [m]', 'Design Wave Height [m]', 'Design Wave Angle [deg]', ...
                'Breaking Water Depth (CEM) [m]', 'Breaking Wave Height (CEM) [m]'};
result_table = table('RowNames', row_names);

for a=1:length(alphas) % For every angle
    parta_table = table();
    for T=1:length(periods) % For every wave period
        [L, ~, ~, ~, alpha, ~, Kr, H] = wave_params(H0, periods(T), d, alphas(a));  % Calculate wavelength, wave angle, and wave height
        H0p = Kr * H0;                                                              % Calculate unrefracted deepwater wave height
        L0 = dispersion(periods(T), d);                                             % Calculate deepwater wavelength
        [db, Hb] = break_params(H0p, H0, L0, periods(T), beta);            % Calculate breaking water depth and breaking wave height
        parta_table = addvars(parta_table, [L; H0p; H; alpha; db; Hb], 'NewVariableNames', col_names{T});
    end
    parta_table = mergevars(parta_table, 1:length(periods), 'NewVariableName', super_col_names{a}, 'MergeAsTable', true);
    result_table = [result_table, parta_table];
end

disp(result_table)

%% Part B
% Depth from Wavelength Observations: Aerial photographs of a coastline display the presence of two wave systems, 
% one with crests 76m apart and wave period of 10 seconds and another with crests 15m apart.
L1 = 76; % m
T1 = 10; % s
L2 = 15; % m

% 1) water depth in which wavelengths were observed:
L0 = 9.81 * T1^2 / (2*pi);
d = L1*atanh(L1/L0) / (2*pi); % water depth [m] - derived from dispersion

% 2) wave period of the shorter wave system:
T2 = sqrt(2*pi*L2 / 9.81*tanh(2*pi/L2*d)); % wave period [s] - derived from dispersion