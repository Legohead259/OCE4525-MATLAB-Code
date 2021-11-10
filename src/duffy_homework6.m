% OCE4525 - Homework 6 - Forces on Vertical Walls
% Braidan Duffy
% 11/12/2021

clear all; close all;

%% Part A
% Water depth at toe of structure ds= 3.0m; offshore slope m = 0.02;
% deep water wave height H0 = 1.50 m;  wave period T =10 sec.;
% deep water wave angle α0 = 0 deg; seawater unit weight ω = 10 kN/m3

ds = 3.0; % m
m = 0.02; % m
beta = atand(m); % deg
H0 = 1.5; % m
T = 10; % s
alpha_0 = 0; % deg
rho_sw = 1025; % kg/m^3
g = 9.81; % m/s/s
W_sw = rho_sw * g; % N/m^3

% 1. Determine the incident wave height at structure
[~, ~, ~, ~, ~, ~, ~, Hi, ~, ~] = wave_params(H0, T, ds, alpha_0);

% 2. For the wave conditions above acting on a solid, vertical wall 
% (with sufficient top elevation to prevent wave overtopping and 100% reflection) at this site, 
% calculate the maximum total wave forces (a) with the crest at the wall and 
% (b)with the trough at the wall, and the total wave moment about the toe of the wall 
% (c)with the crest at the wall and (d)with the trough at the wall, per unit length of wall:

Cr = 1.0; % reflection coefficient

Fs = W_sw * ds^2/2; % Maximum wave force on the wall - N/m
Ms = W_sw * ds^3/6; % Maximum wave moment on the about the toe - N-m/m

p1 = Hi/ds;                 % Dimensionless parameter for SPM Figs
p2 = Hi/(g*T^2);            % Dimensionless parameter for SPM Figs
p3 = 0.915;                 % h0/Hi - Dimensionless parameter from SPM Fig. 7-90;
h0 = p3 * Hi;               % Center of clapotis above MSL - m
yc = ds + h0 + (1+Cr)/2*Hi; % Water depth from the crest of clapotis - m
yt = ds + h0 - (1+Cr)/2*Hi; % Water depth from the crest of clapotis - m
min_height = max([yc, yt]); % Minimum wall height to prevent overtopping - m

p4 = 0.865;             % Fc/Wd^2 - Dimensionless parameter from SPM Fig. 7-91
p5 = -0.375;            % Ft/Wd^2 - Dimensionless parameter from SPM Fig. 7-91
Fc = p4 * W_sw * ds^2;  % Force from clapotis crest - N/m
Ft = p5 * W_sw * ds^2;  % Force from clapotis trough - N/m
Fc_total = Fs + Fc;     % Total force on the vertical wall from the wave crest - N/m
Ft_total = Fs + Ft;     % Total force on the vertical wall from the wave trough - N/m

p6 = 0.65;              % Mc/Wd^2 - Dimensionless parameter from SPM Fig. 7-92
p7 = -0.14;             % Mt/Wd^2 - Dimensionless parameter from SPM Fig. 7-92
Mc = p6 * W_sw * ds^3;  % Moment from clapotis crest - N-m/m
Mt = p7 * W_sw * ds^3;  % Moment from clapotis trough - N-m/m
Mc_total = Ms + Mc;     % Total moment on the vertical wall from the wave crest - N-m/m
Mt_total = Ms + Mc;     % Total moment on the vertical wall from the wave trough - N-m/m

%% Part B
% Given design site data:  water depth at toe of structure ds= 3.0m; offshore slope m = 0.02;
% deep water wave height H0 = 6.00 m;  wave period T =20 sec.;
% deep water wave angle α0 = 0 deg; seawater unit weight ω = 10 kN/m3

clear all;

ds = 3.0; % m
m = 0.02; % m
beta = atand(m); % deg
H0 = 6.00; % m
T = 20; % s
alpha_0 = 0; % deg
rho_sw = 1025; % kg/m^3
g = 9.81; % m/s/s
W_sw = rho_sw * g; % N/m^3

% 1. Determine the incident wave height at structure
[~, ~, ~, ~, ~, ~, ~, Hi, ~, ~] = wave_params(H0, T, ds, alpha_0);
