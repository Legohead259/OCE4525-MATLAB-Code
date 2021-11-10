% OCE4525 - Homework 5 - Forces on Pilings
% Braidan Duffy
% 10/29/2021

clear all; close all;

% Cocoa Beach wants to build a new pier with a maximum depth of 3m below the SWL,  d =3.0m. 
% Calculate the wave forces on the pier pilings using non-linear wave theory. 
% Use the design wave parameters calculated in HW#4 for d = 3m (Hi = 2.97m / T = 20s / breaking wave). 
% Diameter of pilings D = 1.0 m / seawater density = 1025 kg/m3

d = 3; % m
Hi = 2.97; % m
T = 20; % s
D = 1; % m
rho_sw = 1025; % Density - kg/m3
g = 9.81; % Gravitational acceleration - m/s2
mu_sw = 0.00109; % Dynamic viscosity - Ns/m2
[L, C, n, Cg] = wave_params(Hi, T, d); % Calculate wave parameters
[~, ~, u] = wave_particles(Hi, T, d, -d); % Calculate particle velocities

% Question 1
piling_ratio = D/L;
if piling_ratio > 0.05
    disp("Morrison's equation is not valid!")
else
    disp("Morrison's equation is valid!")
end

% Question 2
Re = rho_sw*u*D/mu_sw; % Calculate Reynold's number around the piling

if (Re > 2.52E5 && Re < 5E5) % Calculate Morrison's coefficient
    Cm = 2.5-Re/5E5;
elseif (Re > 5E5)
    Cm = 1.5;
end

if (Re > 2E5 && Re < 5E5) % Calculate drag coefficient
    Cd = 1.2 - (Re-2E5)/6E5;
elseif (Re > 5E5)
    Cd = 0.7;
end

% Question 3 - Note: Use CEM Fig. VI-5-126 thru VI-5-129
Ki = 0.38;
Kd = 0.90;
Si = 0.86;
Sd = 1.15;

% Question 4
Fi = Cm*rho_sw*g*pi*D^2/4*Hi*Ki; % Calculate inertial force on pile
Fd = Cd*1/2*rho_sw*g*D*Hi^2*Kd; % Calculate drag force on pile
Mi = Fi*d*Si; % Calculate maximum inertial moment on pile
Md = Fd*d*Sd; % Calculate maximum drag moment on pile

% Question 5
W = Cm*D/(Cd*Hi);
phi_m = 0.46; % Linearlly interpolated from 0.42 (W=0.5) and 0.5 (W=1.0) gathered from CEM Figs VI-5-133 and 134
alpha_m = 0.41; % Linearlly interpolated from 0.40 (W=0.5) and 0.42 (W=1.0) gathered from CEM Figs VI-5-137 and 138

% Question 6
Fm = phi_m*Cd*rho_sw*g*D*Hi^2; % Calculate combined drag and intertial force on pile
Mm = alpha_m*Cd*rho_sw*g*D*Hi^2*d; % calculate combined drag and inertial moments on pile

% Question 7
KC = u*T/D; % Calculate Keulegan-Carpenter number
if KC > 3
    disp("Lift forces are significant since KC > 3")
else
    disp("Lift forces are insignificant since KC < 3")
end

% Question 8
beta = Re/KC; 
Cl = 0.1; % Lift coefficient as determined by CEM Fig VI-5-141
Fl = Cl*rho_sw*g/2*D*Hi^2*Kd; % Lift force

%% Part B
% Determine the maximum pile diameter (D) for which the Morison Equation is valid for each condition:
% T = [20, 4, 20, 4]
% h = [deepwater, deepwater, 3, 3]

[L1] = dispersion(20, 0); % Calculate first deepwater wavelength
[L2] = dispersion(4, 0); % calculate second deepwater wavelength
[~, L3] = dispersion(20, 3); % Calculate first shallow-water wavelength
[~, L4] = dispersion(4, 3); % calculate second shallow-water wavelength

D1 = L1*0.05; % Calculate maximum pile diameter for T=20s in deepwater
D2 = L2*0.05; % calculate maximum pile diameter for T=4s in deepwater
D3 = L3*0.05; % Calculate maximum pile diameter for T=20s in h=3m
D4 = L4*0.05; % Calculate maximum pile diameter for T=4s in h=3m

%% Part C
% Determine the maximum scour depth (S_m) for the piles and wave conditions in part A, then recalculate the Moment. 
Sm_est = 2*D; % Estimated scour depth
Sm_full = 1.3*(1-exp(-0.03*(KC-6)))*D; % Scour depth as determined by CEM Eq. VI-5-267

ra = Mm / Fm; % Calculate pile moment arm
Mm_prime_est = (ra+Sm_est)*Fm; % Calculate new estimated moment arm, factoring scour
Mm_prime_full = (ra+Sm_full)*Fm; % Calculate new moment arm, factoring scour, from CEM eq.