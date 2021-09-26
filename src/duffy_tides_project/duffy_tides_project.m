% OCE4525 - MP6 - Water Levels and Tides
% Braidan Duffy
% 09/24/2021

clear all; close all;

% Water levels and tides

% Visit www.tidesandcurrents.noaa.gov and select a coastal tide gauge of your choosing. Pick a unique location that others are not already analyzing.

%    30 pts Use the web site to plot the tides for a 1 month time period. Save the plot and copy into word document, on the image you saved, circle and label the spring and neap cycles.
%    10 pts Identify if your site has diurnal, semi-diurnal, or mixed tides
%    10 pts Calculate the difference in max tidal elevation between the spring/neap cycles at that location
%    Extract that historical data to a file and save that file to your computer
%    30 pts Plot the data in the file using MATLAB (be sure to use complete and appropriate labels)
%    20 pts In MATLAB Calculate the mean water level over the duration of your plotted signal; 
%        compare this value to that value given in the information benchmark section of the tide gauge homepage, are they different? why?

% Write up results in MS Word include all figures from web source and plots from MATLAB.  Include MATLAB code as appendix.

% Display NOAA image
figure(1)
imshow(imread('noaanosco-opsobserved-wa.png'))

% API Call parameters
url = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter";
id = "8722670";
name = "West Palm Beach, Atlantic Ocean, FL";
units = "metric";
bdate = "20210825"; % Beginning date - YYYYMMDD
edate = "20210924"; % End date - YYYYMMDD
datum = "mllw"; % Reference datum for measurements
interval = "6"; % Time between measurements - minutes

fn = id + '_' + bdate + '_' + edate + '.csv';

if ~isfile(fn) % Check if the data file has not been created already
    % Make call to NOAA API for station data
    request = webread(url, 'station', id, ...
                        'units', units, ...
                        'begin_date', bdate, ...
                        'end_date', edate, ...
                        'product', 'water_level', ...
                        'time_zone', 'gmt', ...
                        'datum', datum, ...
                        'interval', interval, ...
                        'format', 'json');
    time = extractfield(request.data, 't');                                                 % Get time data
    new_time = datetime(time, 'InputFormat', 'yyyy-MM-dd HH:mm', 'TimeZone', '+00:00');     % Format time data for plotting
    level = extractfield(request.data, 'v');                                                % Get water level data
    new_level = str2double(level);                                                          % Format water level data for plotting
    table = table(new_time, new_level, 'VariableNames', ["Time", "WaterLevel"]);
    

    fid = fopen(fn+'.meta', 'w'); % Open a metadata file
    fprintf(fid, '%s\n', request.metadata.id); % Write the station id
    fprintf(fid, '%s\n', request.metadata.name); % Write the station name
    fprintf(fid, '%s\n', request.metadata.lat); % Write the station latitude
    fprintf(fid, '%s\n', request.metadata.lon); % Write the station longitude
    fclose(fid); % Close the metadata file

    fid = fopen(fn, 'w'); % Open the CSV file
    for i=1:length(time) % For every time value
        fprintf(fid, '%s,', time{i}); % Write the timestamp to file
        fprintf(fid, '%12.3f\n', new_level(i)); % Write the data to file
    end
    fclose(fid); % Close file
else
    opts = delimitedTextImportOptions('NumVariables', 2);
    opts.Delimiter = ',';
    opts.VariableNames = ["Time", "WaterLevel"];
    opts.VariableTypes = ["datetime", "double"];
    opts = setvaropts(opts, 'Time', 'TimeZone', '+00:00', 'InputFormat', 'yyyy-MM-dd HH:mm');
    table = readtable(fn, opts);
end

% State what tidal cycles the site has and Spring/Neap Tides
disp("Site " + id + " at " + name + " has semdiurnal and mixed tidal cycles.")
disp("There is a neap tide around August 31, 2021 with a high tide 0.852m above MLLW")
disp("There is a spring tide around September 8, 2021 with a high tide 1.151m above MLLW")
disp("The difference between neap tide and spring tide between August 25 and September 24 is " + num2str(1.151-0.852) + "m above MLLW")
disp(" ")

% Calculate mean tide level and compare
MTL_ref = 0.463; % m re: MLLW
MTL = mean(table.WaterLevel); % m re: MLLW
disp("The mean tide level for August 24 to September 25 is " + num2str(MTL, 3) + "m above MLLW")
disp("The mean tide level datum is 0.463m above MLLW")
disp("The difference between the measured MTL and datum MTL is " + num2str(MTL-MTL_ref) + "m")
disp("This suggests that the water level in the area is higher than normal, possibly due to warmer, less dense water")

% Plot water elevation for the site from NOAA data
figure(2)
plot(table.Time, table.WaterLevel)
xlabel("Time")
ylabel("Water Level (re: " +  datum + ") [m]")
title("Water level at " + name + " (Station ID: " + id + ")")
xline(datetime(2021, 08, 31, 'TimeZone', '+00:00'), '--', 'Neap Tide', 'Color', 'r') % Draw line at neap tide
xline(datetime(2021, 09, 07, 'TimeZone', '+00:00'), '--', 'Spring Tide', 'Color', 'r') % Draw line at spring tide
yline(MTL_ref, '--', 'MTL Datum', 'Color', 'm') % Draw line at reference MTL
yline(MTL, '--', 'Measured MTL', 'Color', 'm') % Draw line at measured MTL