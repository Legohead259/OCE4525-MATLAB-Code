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

L_arr = zeros(length(periods), length(alphas));         % Instantiate matrix of wavelengths
H0p_arr = zeros(length(periods), length(alphas));       % Instantiate matrix of unrefracted deepwater wave heights
Hd_arr = zeros(length(periods), length(alphas));        % Instantiate matrix of wave heights at specified water depth
alpha_arr = zeros(length(periods), length(alphas));     % Instantiate matrix of wave angles
db_arr = zeros(length(periods), length(alphas));        % Instantiate matrix of breaking water depths (CEM method)
Hb_arr = zeros(length(periods), length(alphas));        % Instantiate matrix of breaking wave heights (CEM method)

for a=1:length(alphas) % For every angle
    for T=1:length(periods) % For every wave period
        [L_arr(T,a), ~, ~, ~, alpha_arr(T,a), ~, Kr, Hd_arr(T,a)] = wave_params(H0, periods(T), d, alphas(a));  % Calculate wavelength, wave angle, and wave height
        H0p_arr(T,a) = Kr * H0;                                                                                 % Calculate unrefracted deepwater wave height
        L0 = dispersion(periods(T), d);                                                                         % Calculate deepwater wavelength
        [db_arr(T,a), Hb_arr(T,a)] = break_params(H0p_arr(T,a), H0, L0, periods(T), beta);                      % Calculate breaking water depth and breaking wave height
    end
end

outputTable = {length(alphas)};
for i=1:length(alphas) % For every angle, create a corresponding table
    outputTable{i} = table(L_arr(:, i), H0p_arr(:, i), Hd_arr(:, i), alpha_arr(:, i), db_arr(:, i), Hb_arr(:, i));
    outputTable{i}.Properties.VariableNames = {'Wavelength [m]', 'Unrefracted H0 [m]', 'Wave height at d [m]', ...
                                                'Wave angle at d [deg]', 'Breaking water depth [m]', 'Breaking wave height [m]'};
    % Display
    disp(outputTable{i})
end

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