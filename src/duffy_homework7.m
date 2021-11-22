% OCE4525 - Homework 7 - Rubble Mound Breakwater Design
% Braidan Duffy
% Collaborated with Gabbriella S. and Clay K.
% 11/22/2021

clear all; close all;

%% Introduction
% Site Conditions: Incident breaking wave height = 5.8 m, wave period = 20 sec., 
% water depth = 6m, seabed slope = 0.02, seawater density = 1,025 kg/m3 .
H0 = 5.8; % m
T = 20; % sec
h = 6; % m
m = 0.02;
beta = atand(m); % deg
rho_sw = 1025; % kg/m3
g = 9.81; % m/s/s
L0 = g*T^2/(2*pi); % m

%% Part A
% 1. Rubble Mound Quarrystone Structure Design:  Calculate the armor unit sizes and design 
% details for a rubble mound structure at this site with no overtopping allowed, 
% on a 1V:2H structure trunk slope, using 2 layers of randomly placed rough angular 
% quarrystone granite rock with a mass density = 2,650 kg/m3:

parta_ans = table('RowNames', ["Mean unit mass (M_50) [kg]", "Mean unit diameter (D_n50) [m]", ...
                   "Armor layer thickness (r) [m]", "Armor unit density (N/A) [units/m^2]", ...
                   "Average wave runup (R_ave) [m]"]);

rho_s = 2650;
delta = (rho_s / rho_sw) - 1;
a = .013;
b = 22;
Hs = H0;
cota = 2;
Kd1 = 2.0;
M_n50 = (rho_s * Hs ^ 3) / (Kd1 * (delta ^ 3) * cota);
D_n50 = (M_n50 / rho_s)^(1/3);
n1 = 2;
kdell = 1;
r = n1 * kdell * D_n50;
P1 = 37;
gamma_1 = rho_s * g;
Na_A = n1 * kdell * (rho_s/gamma_1)^(2/3) * (1 - P1/100);
tanb = .03;
delta = (rho_s / rho_sw) - 1;
R1 = Hs / (delta*D_n50);
Ns1 = R1 * (Hs / L0)^(-1/3);
R_ave = Hs * tanb / sqrt(Hs / L0);

answers = [M_n50; D_n50; r; Na_A; R_ave];
parta_ans = addvars(parta_ans, answers, 'NewVariableNames', 'No Allowable Damage (Quarrystone)');

Hs25 = H0 / 1.37;
Kd = 31;
cota2 = 2;
M_n50 = (rho_s * Hs25 ^ 3) / (Kd * (delta ^ 3) * cota2);    
D_n50 = (M_n50 / rho_s)^(1/3);
r = n1 * kdell * D_n50;
W_50 = M_n50 * g;
Na_A = n1 * kdell * (gamma_1/W_50)^(2/3) * (1 - P1/100);
R_ave = Hs * tanb / sqrt(Hs25 / L0);

answers = [M_n50; D_n50; r; Na_A; R_ave];
parta_ans = addvars(parta_ans, answers, 'NewVariableNames', '25% Allowable Damage (Quarrystone)');

%% Part A Section 2
% Rubble Mound Concrete Armor Unit Structure Design:  
% Calculate the armor unit sizes and design details for a rubble mound structure at this site 
% with no overtopping allowed, on a 1V:2H structure trunk slope, using 2 layers of randomly placed 
% concrete tetrapods with a mass density = 2,400 kg/m3:

rho_s = 2400;
delta = (rho_s / rho_sw)-1;        
Kd = 7;
M_n50 = (rho_s * H0 ^ 3) / (Kd * (delta ^ 3) * cota);
D_n50 = (M_n50 / rho_s)^(1/3);
kdell = 1.04;
r = n1 * kdell * D_n50;
W_50 = M_n50 * g;
gamma = rho_s * g;
P3 = 50;
Na_A = n1 * kdell * (gamma/W_50)^(2/3) * (1 - P3/100);
R_ave = Hs * tanb / sqrt(Hs / L0);

answers = [M_n50; D_n50; r; Na_A; R_ave];
parta_ans = addvars(parta_ans, answers, 'NewVariableNames', 'No Allowable Damage (Tetrapods)');

Hs_25 = Hs / 1.32;
Kd = 7;
M_n50 = (rho_s * Hs_25 ^ 3) / (Kd * (delta ^ 3) * cota);    
D_n50 = (M_n50 / rho_s)^(1/3);
r = n1 * kdell * D_n50;
W_50 = M_n50 * g;
Na_A = n1 * kdell * (gamma/W_50)^(2/3) * (1 - P3/100);
R_ave = Hs * tanb / sqrt(Hs_25 / L0);

answers = [M_n50; D_n50; r; Na_A; R_ave];
parta_ans = addvars(parta_ans, answers, 'NewVariableNames', '25% Allowable Damage (Tetrapods)');

%% Part B - Van der Meer Formula Calculations
% 1. Rubble Mound Quarrystone Structure Design:  Calculate the armor unit sizes and 
% design details for a rubble mound structure at this site with no overtopping allowed, 
% on a 1V:2H structure trunk slope, using 2 layers of randomly placed rough angular 
% quarrystone granite rock with a mass density = 2,650 kg/m3:

Nz = 7500;                      % Number of waves designed for 
P = 0.4;                        % Porosity
n = 2;                          % Number of armor stone layers
kdell = 1;                      % Layer coefficient from CEM Table VI-5-51
rho_s = 2650;                   % Quarrystone density - kg/m
w_a = rho_s * g;                % Specific gravity of armor stones
delta = (rho_s / rho_sw) - 1;
Hs = H0;                        % Design significant wave height is incident H - m
m_BW = 1/2;                     % Slope of the BW face
alpha = atand(m_BW);            % Angle of the breakwater face - deg

partb_ans = table('rownames', ["Mean unit mass (M_50) [kg]", "Mean unit diameter (D_n50) [m]", ...
                   "Armor layer thickness (r) [m]", "Armor unit density (N/A) [units/m^2]", ...
                   "2% wave runup (R_2p) [m]"]);

% No Allowable Damage (D=0%)
Kd = 4;                                         % Stability coefficient from SPM Table 7-7
S = 2;                                          % Damage level from CEM Table VI-5-21

xi_m = m_BW / sqrt(Hs/L0);                      % Iribarren number
xi_mc = ((6.2*P^0.31)/sqrt(m_BW))^(1/(P+1/2));  % Critical Iribarren number
if xi_m < xi_mc % Comapare Iribarren and critical Iribarren numbers to determine mean rock diameter
    D_n50 = Hs / (delta*6.2*S^0.2*P^0.18*Nz^-0.1*xi_m^-0.5);
else
    D_n50 = Hs / (1*S^0.2*P^-0.13*Nz^-0.1*cot(alpha)^0.5*xi_m^P);
end

M_n50 = rho_s * 2.53^3; % Approximate mass of rocks by assuming they're cubes

r = n * kdell * (M_n50/rho_s)^(1/3); % Armor layer thickness - m - CEM Eq. VI-5-117

Na_A = n * kdell * (1-P/100) * (rho_s/M_n50)^(2/3); % Placing density - units/m2 - CEM Eq. VI-5-118

% Runup Coefficients - CEM Table VI-5-5
A = 0.96;
B = 1.17;
C = 0.46;
D = 1.97;

if xi_m > 1 && xi_m <= 1.5
    R_2p = Hs * A * xi_m;
elseif xi_m > 1.5 && xi_m <= (D/B)^(1/C)
    R_2p = Hs * B * xi_m^C ;
elseif xi_m > (D/B)^(1/C) && xi_m < 7.5
    R_2p = Hs * D;
end

answers = [M_n50; D_n50; r; Na_A; R_2p];
partb_ans = addvars(partb_ans, answers, 'NewVariableNames', 'No Allowable Damage (Quarrystone)');

% 25% Allowable Damage
Kd = 10;                                         % Stability coefficient from SPM Table 7-7
S = 8;                                          % Damage level from CEM Table VI-5-21

xi_m = m_BW / sqrt(Hs/L0);                      % Iribarren number
xi_mc = ((6.2*P^0.31)/sqrt(m_BW))^(1/(P+1/2));  % Critical Iribarren number
if xi_m < xi_mc % Comapare Iribarren and critical Iribarren numbers to determine mean rock diameter
    D_n50 = Hs / (delta*6.2*S^0.2*P^0.18*Nz^-0.1*xi_m^-0.5);
else
    D_n50 = Hs / (1*S^0.2*P^-0.13*Nz^-0.1*cot(alpha)^0.5*xi_m^P);
end

M_n50 = rho_s * 1.91^3; % Approximate mass of rocks by assuming they're cubes

r = n * kdell * (M_n50/rho_s)^(1/3); % Armor layer thickness - m - CEM Eq. VI-5-117

Na_A = n * kdell * (1-P/100) * (rho_s/M_n50)^(2/3); % Placing density - units/m2 - CEM Eq. VI-5-118

% Runup Coefficients - CEM Table VI-5-5
A = 0.96;
B = 1.17;
C = 0.46;
D = 1.97;

if xi_m > 1 && xi_m <= 1.5
    R_2p = Hs * A * xi_m;
elseif xi_m > 1.5 && xi_m <= (D/B)^(1/C)
    R_2p = Hs * B * xi_m^C ;
elseif xi_m > (D/B)^(1/C) && xi_m < 7.5
    R_2p = Hs * D;
end

answers = [M_n50; D_n50; r; Na_A; R_2p];
partb_ans = addvars(partb_ans, answers, 'NewVariableName', '25% Allowable Damage (Quarrystone)');

%%
% 2. Rubble Mound Concrete Armor Unit Structure Design: Calculate the armor 
% unit sizes and design details for a rubble mound structure at this site 
% with no overtopping allowed, on a 1V:2H structure trunk slope, using 2 
% layers of randomly placed concrete tetrapods with a mass density =
% 2,400 kg/m3. Non-depth limited waves.

Nz = 7500;                      % Number of waves designed for 
P = 0.4;                        % Porosity
n = 2;                          % Number of armor stone layers
kdell = 1;                      % Layer coefficient from CEM Table VI-5-51
rho_s = 2400;                   % Quarrystone density - kg/m
w_a = rho_s * g;                % Specific gravity of armor stones
delta = (rho_s / rho_sw) - 1;
Hs = H0;                        % Design significant wave height is incident H - m
m_BW = 1/2;                     % Slope of the BW face
alpha = atand(m_BW);            % Angle of the breakwater face - deg

% No Allowable Damage (D=0%)
Kd = 8.3;                                         % Stability coefficient from SPM Table 7-7
S = 2;                                          % Damage level from CEM Table VI-5-21

xi_m = m_BW / sqrt(Hs/L0);                      % Iribarren number
xi_mc = ((6.2*P^0.31)/sqrt(m_BW))^(1/(P+1/2));  % Critical Iribarren number
if xi_m < xi_mc % Comapare Iribarren and critical Iribarren numbers to determine mean rock diameter
    D_n50 = Hs / (delta*6.2*S^0.2*P^0.18*Nz^-0.1*xi_m^-0.5);
else
    D_n50 = Hs / (1*S^0.2*P^-0.13*Nz^-0.1*cot(alpha)^0.5*xi_m^P);
end

M_n50 = rho_s * 2.53^3; % Approximate mass of rocks by assuming they're cubes

r = n * kdell * (M_n50/rho_s)^(1/3); % Armor layer thickness - m - CEM Eq. VI-5-117

Na_A = n * kdell * (1-P/100) * (rho_s/M_n50)^(2/3); % Placing density - units/m2 - CEM Eq. VI-5-118

% Runup Coefficients - CEM Table VI-5-5
A = 0.96;
B = 1.17;
C = 0.46;
D = 1.97;

if xi_m > 1 && xi_m <= 1.5
    R_2p = Hs * A * xi_m;
elseif xi_m > 1.5 && xi_m <= (D/B)^(1/C)
    R_2p = Hs * B * xi_m^C ;
elseif xi_m > (D/B)^(1/C) && xi_m < 7.5
    R_2p = Hs * D;
end

answers = [M_n50; D_n50; r; Na_A; R_2p];
partb_ans = addvars(partb_ans, answers, 'NewVariableNames', 'No Allowable Damage (Concrete)');

% 25% Allowable Damage
Kd = 19.2;                                         % Stability coefficient from SPM Table 7-7
S = 8;                                          % Damage level from CEM Table VI-5-21

xi_m = m_BW / sqrt(Hs/L0);                      % Iribarren number
xi_mc = ((6.2*P^0.31)/sqrt(m_BW))^(1/(P+1/2));  % Critical Iribarren number
if xi_m < xi_mc % Comapare Iribarren and critical Iribarren numbers to determine mean rock diameter
    D_n50 = Hs / (delta*6.2*S^0.2*P^0.18*Nz^-0.1*xi_m^-0.5);
else
    D_n50 = Hs / (1*S^0.2*P^-0.13*Nz^-0.1*cot(alpha)^0.5*xi_m^P);
end

M_n50 = rho_s * 1.91^3; % Approximate mass of rocks by assuming they're cubes

r = n * kdell * (M_n50/rho_s)^(1/3); % Armor layer thickness - m - CEM Eq. VI-5-117

Na_A = n * kdell * (1-P/100) * (rho_s/M_n50)^(2/3); % Placing density - units/m2 - CEM Eq. VI-5-118

% Runup Coefficients - CEM Table VI-5-5
A = 0.96;
B = 1.17;
C = 0.46;
D = 1.97;

if xi_m > 1 && xi_m <= 1.5
    R_2p = Hs * A * xi_m;
elseif xi_m > 1.5 && xi_m <= (D/B)^(1/C)
    R_2p = Hs * B * xi_m^C ;
elseif xi_m > (D/B)^(1/C) && xi_m < 7.5
    R_2p = Hs * D;
end

answers = [M_n50; D_n50; r; Na_A; R_2p];
partb_ans = addvars(partb_ans, answers, 'NewVariableName', '25% Allowable Damage (Concrete)');