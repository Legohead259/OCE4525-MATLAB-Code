% OCE4525 - Wave Tool
% Braidan Duffy
% 08/31/2021

close all; clear all;

% Constants
g = 9.81;       % m/s^2
rho_sw = 1025;  % kg/m^3

% User Parameters
d1 = 2;        % Depth at position 1 - m
H1 = 1.5;         % Wave height at d_1 - m
theta1 = 0;    % Wave angle at d_1 - deg rel to shore
T = 12;         % Wave period - s
m = 1/50;       % Bottom slope
d2 = 5;         % Depth at position 2 - m
z2 = -2.5;      % Depth below SWL at position 2 - m

% Wave characteristics at d1
[L1, C1, n1, Cg1, ~, Ks1, ~, H1, E1, F1] = wave_params(H1, T, d1); % Calculate the wave parameters at d_1

% Deepwater characteristics
[L0, ~, k0, k1] = dispersion(T, d1);    % Get deepwater wavelength and wave numbers at d_0 and d_1
theta0 = asind(k1*sind(theta1)/k0);     % Get the deep water shore angle from Snell's law
Kr1 = sqrt(cosd(theta0)/cosd(theta1));  % Get the refraction coefficient
H0 = H1/Ks1/Kr1;                        % Determine the deepwater wave height (derived from energy flux)
C0 = L0 / T;                            % Get the deepwater celerity
Cg0 = 0.5 * C0;                         % Get the deepwater group celerity
H0p = H1/Ks1;                           % Get the deepwater waveheight without factoring in refraction

% Breaking depth characteristics
[Hb, db] = break_params(H0p, H0, L0, T, atand(m));

if d2 > db % Wave does not break before d2, calculate wave characteristics and particle movement at d2
    [L2, C2, ~, Cg2, theta2, Ks2, Kr2, H2, E2, F2] = wave_params(H0, T, d2, theta0);    % Calculate the wave parameters at d_1
    E2 = E2 * L2;                                                                     % Calculate energy per unit crest width [J/m]
    H02 = H2 / Ks2 / Kr2;
    H02p = H02 * Kr2;
    [Sx, Sz, u, w, ax, az] = wave_particles(H2, T, d2, z2);                         % Calculate particle movements
else % Wave breaks before d2
    disp("Wave breaks before d2. Wave characteristics at d2 not calculated")
end
    