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
m = 0.02;
beta = atand(m); % deg
H0 = 6.00; % m
T = 20; % s
alpha_0 = 0; % deg
rho_sw = 1025; % kg/m^3
g = 9.81; % m/s/s
W_sw = rho_sw * g; % N/m^3

% 1. Determine the incident wave height at structure
[~, ~, ~, ~, ~, ~, ~, Hi, ~, ~] = wave_params(H0, T, ds, alpha_0);

p1 = ds / (g*T^2);  % Dimensionles parameter for SPM Fig. 7-4
p2 = .99;           % Hb/ds - From SPM Fig. 7-4
Hb = ds * p2;       % Determine breaking wave height from SPM Fig. 7-4

[L0, L] = dispersion(T, ds);
Ph = rho_sw * g * (ds+(Hb/2));
hc = .78 * Hb;
rs = rho_sw * g * .5 * (ds + Hb/2)^2;
ms = (rho_sw * g * (ds + Hb/2)^3)/6;
pm_fig = 101*rho_sw*g*(Hb/L) *(ds/(ds+L))*(ds+L+ds);

f1 = ds/L0;
f2 = ds/L;
Ld = ds/f2;
db  = ds+Ld;
f3 = db/L0;
f4 = db/Ld;

db = Ld + ds;
f5 = db/L0;
f6 = db/L0;
f7 = db/f6;
pm_eqn = 3.5*rho_sw*g*Hb;
rm_eqn = 3.5*rho_sw*g*Hb^2 / 3;
rm_fig = pm_fig*Hb/3;
Mm_fig = rm_fig * ds;
Mm_eqn = pm_eqn *Hb*ds/3;

hc = Hb/2;                  % Breaking wave crest height
Ps = W_sw*(ds + hc);           % Static pressure on the wall - Pa
Rs = (W_sw*(ds + Hb/2)^2)/2;   % Static force on the wall - N/m
Ms = (W_sw*(ds + Hb/2)^3)/6;   % Static moment on the wall - N-m/m

db = 5.15;                              % Breaking depth of the wave - m
[~, Lb] = dispersion(T, db);            % Breaking wavelength - m
Pm = 101*(W_sw)*(Hb/Ld)*(ds/db)*(db+ds);  % Dynamic pressure on the wall - Pa
Rm = (Pm*Hb)/3;                         % Dynamic force on the wall - N/m
Mm = (Pm*Hb*ds)/3;                      % Dynamic moment on the wall - N-m/m

Pt = Ps + Pm; % Total pressure on the wall - Pa
Rt = Rm + Rs; % Total force on the wall - N/m
Mt = Mm + Ms; % Total moment on the wall - N-m/m

% fig-7-100 data
p3 = ds/(g*T^2); % Dimensionless x-parameter for SPM Fig. 7-100
p4 = 3; % Pm/(W*Hb) from SPM Fig. 7-100
Pm_fig = W_sw*Hb*p4; % Pressure from SPM Fig. 7-100 - Pa
Rm_fig = (1/3)*W_sw*Hb^2*3; % Force from SPM Fig. 7-100 - N/m
Mm_fig = (Pm_fig*Hb*ds)/3; % Moment from SPM Fig. 7-100 - N-m/m
Pt_fig = Ps + Pm_fig; % Total pressure from SPM Fig. 7-100 - Pa
Rt_fig = Rm_fig + Rs; % Total force from SPM Fig. 7-100 - N/m
Mt_fig = Mm_fig + Ms; % Total moment on the wall from SPM Fig. 7-100 - N-m/m

% Determine min height of wall (assuming perfect reflection)
Cr=1;
y_c = ds + ((1+ Cr)/2)*Hb; 

%% Part C
% Concrete vertical wall Caisson resting on the seabed; width of the caisson B = 3.0m; 
% water depth at toe of structure ds= 3.0m; offshore slope m = 0.02; 
% deep water wave height H0 = 1.50 m; wave period T =10 sec.; 
% deep water wave angle beta_0 = 0 deg; seawater unit weight ω = 10 kN/m3

clear all;

B = 3.0; % m
ds = 2.0; % m
m = 0.02;
beta = atand(m); % deg
H0 = 1.5; % m
T = 10; % s
beta_0 = 0; % deg
rho_sw = 1025; % kg/m^3
rho_c = 2400; % Density of concrete - kg/m^3
g = 9.81; % m/s/s
W_sw = rho_sw * g; % N/m^3

% 1. Determine the incident wave height at structure
[~, ~, ~, ~, ~, ~, ~, Hi, ~, ~] = wave_params(H0, T, ds, beta_0);
[~, ~, ~, ~, ~, kh] = dispersion(T, ds);
% TODO: Determine minimum wall height to prevent over-topping
hw = 5.78; % Minimum wall height to prevent over-topping - m

% For the wave conditions above acting on a vertical wall caisson 
% (with sufficient top elevation to prevent wave overtopping and 100% reflection) 
% at this site, calculate the maximum total wave forces and the total wave moments 
% about the toe of the wall.  Use the Goda Method outlined in the CEM to solve problem, 
% and compare your answer to Part A above.

% Determine pressure due to the wave on the wall
alpha_1 = 0.6 + 0.5*(2*kh/sinh(2*kh))^2; % Wave period effect
db = ds + m*5*Hi; % Water depth 5Hs seaward of the structure face - m
a1 = (db-ds)/(3*db) * (Hi/ds)^2; % Factor 1 for alpha_2
a2 = 2*ds/Hi; % Factor 2 for alpha_2
alpha_2 = min([a1 a2]); % Increase in pressure due to mound
p1 = 0.5 * (1+cosd(beta_0)) * (alpha_1 + alpha_2*cosd(beta_0)^2) * W_sw * Hi; % Pressure factor 1 - Pa

% Determine pressure at the top of the wall
Hc = Hi/2;
eta_star = 0.75*(1+cos(beta_0)) * Hi; % Overtopping wave height - m
if eta_star > Hc
    p2 = (1-Hc/eta_star); % Pressure factor 2 - Pa
else
    p2 = 0; % Pressure factor 2 - Pa
end

% Determine pressure at the tow of the wall
alpha_3 = 1 - (hw-Hc)/ds * (1-1/cosh(kh)); % Pressure distribution across face
p3 = p1*alpha_3; % Pressure factor 3 - Pa

% Determine pressure on the bottom of the wall
pu = 0.5 * (1+cosd(beta_0)) * alpha_1 * alpha_3 * W_sw * Hi; % Pressure factor u - Pa

% Determine Forces on the wall
% Uncertainty factors and bias of horizontal wave induced force, uplift force, horizontal moment, and uplift moment
UFH = 0.90;
UFU = 0.77;
UMH = 0.81;
UMU = 0.72;

FH = UFH * (0.5*(p1+p2)*Hc + 0.5*(p1+p3)*ds);   % Horizontal force on the wall - N/m
FU = UFU * 0.5*pu * B;                          % Uplift force on the wall - N/m
FG = rho_c*g*B*hw - rho_sw*g*B*ds;              % Some force on the wall - N/m

MH = FH * ds/2; % Horizontal moment on the wall - N-m/m
MU = FU * Hi; % Uplift moment on the wall - N-m/m
MG = FG * db; % Some moment on the wall - N-m/m