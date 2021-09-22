% OCE4525 - Homework 3: Diffraction
% Braidan Duffy
% 09/22/21

% A long straight breakwater runs from south to north into an area of the
% ocean where the depth landward (west) of the breakwater is constant, d = 10m.
% Design waves from NE, E, and SE directions with periods T1 = 8 seconds 
% and T2 =16 seconds and incident wave height,  Hi = 3.0m,  are arriving 
% at the breakwater at incident angles of 45, 90 and 120 degrees azimuth from north.
% A berth for ships is planned for a point that is located 500 meters on a 
% radius of 45 degrees inside of the tip of the breakwater (SW of breakwater 
% tip at an azimuth of  225 degrees).

% Find the wave heights at the ship berth for the six different wave conditions, 
% and the deep-water wave conditions that would have generated the incident wave conditions, 
% assuming straight and parallel offshore depth contours:

% Constants
d = 10; % Landward water depth - m
T = [8,16]; % Design wave periods - s
Hi = 3.0; % Incident wave height - m
alphai = [45, 0, -30]; % Incident wave angles (azimuth from east) - deg
r = 500; % Radius from breakwater to berth - m
theta = 135; % Angle from breakwater to berth (azimuth from east) - deg

% Values from SPM plates
H0p_SPM = Hi ./ [0.9, 0.95]; % Values gathered from SPM figure C-1
Kprime = [0.55, 0.55; 0.115, 0.175; 0.0825, 0.12]; % Values gathered from SPM figures [angle, period]

% Results
super_col_names = {'α=45°', 'α=90°', 'α=120°'};
col_names = {'T=8s', 'T=16s'};
row_names = {'Incident Wavelength [m]', 'Unrefracted H0 (SPM) [m]', 'Unrefracted H0 (LWT) [m]', ...
                'Wave Height at Ship Berth [m]', 'Deepwater Wave Height [m]', 'Deepwater wave angle [deg]'};
result_table = table('RowNames', row_names);

% Calculations
for a=1:length(alphai)
    temp_table = table();
    for t=1:length(T)
        [L0, Li, k0, ki, ~, khi] = dispersion(T(t), d); % Calculate Li and other parameters
        C = Li/T(t);                                    % Calculate wave celerity
        n = (1/2*(1+2*khi/sinh(2*khi)));                % Calculate design group wave modulation
        Cgi = n*C;                                      % Calculate design group wave celerity
        Ks = sqrt(L0/(2*T(t)*Cgi));                     % Calculate shoaling coefficient
        H0p_LWT = Hi/Ks;                                % Calculate unrefracted DW wave height using LWT
        Hd = Hi * Kprime(a,t);                          % Calculate the wave height at the ship berth
        alpha0 = real(asind(ki*sind(alphai(a))/k0));    % Calculate the DW wave angle (real value is taken since alpha0 can exceed α<90°)
        Kr = sqrt(cosd(alpha0)/cosd(alphai(a)));        % Calculate refraction coefficient
        H0 = Hi/Ks/Kr;                                  % Calculate DW wave height
        temp_table = addvars(temp_table, [Li; H0p_SPM(t); H0p_LWT; Hd; H0; alpha0], ...
                                'NewVariableNames', ...
                                col_names{t});          % Add variables into table
    end
    temp_table = mergevars(temp_table, 1:length(T), 'NewVariableName', super_col_names{a}, 'MergeAsTable', true); % Create subtables
    result_table = [result_table, temp_table]; % Concatenate subtables into main result table
end

% Display
disp(result_table)