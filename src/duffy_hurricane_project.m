% OCE4525 - Mini Project 5 - Hurricane Wave Generation
% Braidan Duffy
% Collaborated with: Clay Kenney, Gabriella Santiago
% 09/15/2021

% Part A): Write a MATLAB code to generate the vortex wind fields add the forward velocity 
% vector to the wind field and plot the wind fields at the two indicated times below 
% (48 and 24 hours before landfall).  For both of the hindcast times indicated below, and using 
% the associated NOAA data, construct wind fields for Hurricane Katrina at the instrument height 
% and plot the results using MATLAB. Use both the SPM (Equation 3-57 page 3-82 of SPM Volume 1) 
% and the CEM (Holland Wind model Eq. II-2-18 in CEM) equations and compare your results.

% Part B): Then estimate the wave parameters (H, T) for each case based on the wave growth formulas 
% in the CEM. Using the data from Part 1), estimate the deep water wave parameters using CEM.  
% Assume fetch is equal to the radius of maximum wind and duration is 10 hours.  
% Is the condition fetch limited or duration limited?

% EC: Estimate the water depth based on the position of the storm and shoal the deep water wave 
% calculated in part 2 to the depth estimated here.

clear all; close all;

% Constants
dates = ["12:00 UTC 27 Aug 2005", "12:00 UTC 28 Aug 2005"];
Pc = [94200, 90800];                % Central pressures - Pa
Vf = [3.09, 5.14];                  % Forward velocities - m/s
theta = [275, 300];                 % Forward direction - deg
lat = [24.5, 25.7];                 % Latitude - deg
Pn = 101400;                        % Ambient pressure - Pa
R = 40000;                          % Radius to maximum wind - m
T_air = 295.15;                     % Air temperature - k
rho_air = Pn / (287.058 * T_air);   % Air density - kg/m^3
B = 1.5;                            % Define arbitrary wind speed distribution peakedness parameter for CEM method (between 1.4 and 2.4)
A = R^B;                            % Define scaling parameter for CEM method
X = R;                              % Fetch - m
t = 10*3600;                        % Wind duration - sec
r_arr = 0:1000:500000;              % Create an array of radiuses from 0 to 400 km in 1 km increments - m
phi_arr = 0:1:360;                  % Create an array of angles around a circle from [0, 360) - deg
g = 9.81;                           % Gravitational acceleration - m/s/s

partb_table = table('RowNames', ["Fetch-limited Wave Period [s]", "Fetch-limited Wave Height [m]",...
                                 "Duration Equivalent Fetch [m]", "Duration-limited Wave Height [m]",...
                                 "Duration-limited Wave Period [s]", "Fetch or Duration Limited?"]);

% Calculations
for i=1:length(Pc) % For each hurricane
    f = 2 * 7.292E-5 * sin(lat(i)); % Determine Coriolis parameter
    for j = 1:length(r_arr) % For each radius from the center of the hurricane
        for k = 1:length(phi_arr) % For each angle around the center of the hurricane  
            x(j,k) = r_arr(j) * cosd(phi_arr(k));
            y(j,k) = r_arr(j) * sind(phi_arr(k));
            
            % ---SPM Method---
            spm_P(i,j,k) = exp(-R/r_arr(j))*(Pn-Pc(i)) + Pc(i);
            spm_Vgr(i,j,k) = -r_arr(j)*f/2 + sqrt((r_arr(j)*f/2)^2 + (Pn-Pc(i)/rho_air) * R/r_arr(j)*exp(-R/r_arr(j)));
            spm_V(i,j,k) = 0.865 * spm_Vgr(i,j,k) + 0.5*Vf(i)*cosd(phi_arr(k) + theta(i));
            
            % ---CEM Method---
            cem_P(i,j,k) = Pc(i) + (Pn-Pc(i))*exp(-A/r_arr(j)^B);
            cem_Vgr(i,j,k) = sqrt((r_arr(j)*f/2)^2 + A*B*(Pn-Pc(i))/(rho_air*r_arr(j)^B) * exp(-A/r_arr(j)^B)) - (r_arr(j)*f/2);
            cem_V(i,j,k) = cem_Vgr(i,j,k) + Vf(i)*cosd(theta(i)+phi_arr(k));
        end
    end
    
    % ---Wave Height and Periods---
    cem_Vmax = sqrt((R*f/2)^2 + A*B*(Pn-Pc(i))/(rho_air*R^B) * exp(-A/R^B)) - (R*f/2);
    Z = 1;
    U10 = cem_Vmax * (10/Z)^(1/7);
    Cd = 0.001 * (1.1+0.035*U10);
    Ustar = U10 * sqrt(Cd);

    X_T = Ustar/g * 0.651 * (g*X/Ustar^2)^(1/3);
    X_H = Ustar^2/g * (4.13E-2)*sqrt(g*X/Ustar^2);

    Xt = 5.23E-3 * (g*t/Ustar)^(3/2) * Ustar^2/g;
    t_T = Ustar/g * 0.651 * (g*Xt/Ustar^2)^(1/3);
    t_H = Ustar^2/g * (4.13E-2)*sqrt(g*Xt/Ustar^2);
    
    partb_table = addvars(partb_table, [X_T; X_H; Xt; t_T; t_H; "Fetch-Limited"], 'NewVariableNames', dates(i));
end

disp(partb_table)

% ---Plotting---
plotargs = {'RadialRange',[0 max(r_arr)/1000], 'RadLabels',5, 'RadLabelLocation',{180,'top'}, 'plottype','contour', 'PolarDirection','cw'}; % Create standard set of plot arguments 

% Plot wind field (SPM)
figure
spm_V_plot(:,:) = spm_V(1,:,:);
cl = round(min(spm_V_plot)-1):3:round(max(spm_V_plot)+1);
[~,~,~,c,h] = polarplot3d(spm_V_plot,'contourlines',cl,plotargs{:});
title("Hurricane Katrina Wind Field - " + dates(i) + " (SPM Method)")
clabel(c,h)
legend("48 Hours before landfall")
zlabel("Wind Speed (m/s)")

% Plot pressure field (SPM)
figure
spm_P_plot(:,:) = spm_P(1,:,:);
[~,~,~,c,h] = polarplot3d(spm_P_plot, plotargs{:});
title("Hurricane Katrina Pressure Field - " + dates(i) + " (SPM Method)")
clabel(c,h)
legend("48 Hours before landfall")
zlabel("Atmospheric Pressure (Pa)")

% Plot wind field (CEM)
figure
cem_V_plot(:,:) = cem_V(1,:,:);
cl = round(min(cem_V_plot)-1):5:round(max(cem_V_plot)+1);
[~,~,~,c, h] = polarplot3d(cem_V_plot,'contourlines',cl,plotargs{:});
title("Hurricane Katrina Wind Field - " + dates(i) + " (CEM Method)")
clabel(c,h)
legend("24 Hours before landfall")
zlabel("Wind Speed (m/s)")

% Plot pressure field (CEM)
figure
cem_P_plot(:,:) = cem_P(1,:,:);
[~,~,~,c,h] = polarplot3d(cem_P_plot, plotargs{:});
title("Hurricane Katrina Pressure Field - " + dates(i) + " (CEM Method)")
clabel(c,h)
legend("24 Hours before landfall")
zlabel("Atmospheric Pressure (Pa)")