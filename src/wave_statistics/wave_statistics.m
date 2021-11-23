% OCE4525 - Design Wave Calculations
% Braidan Duffy
% October 13, 2021
% Updated: November 18, 2021

close all; clear all;

%% Download data
% Historical data
start_year = 1988;
end_year = 2020;
N_years = end_year - start_year + 1;
station_id = "41009"; % NOAA NDBC station number

for i=start_year:end_year % For every historical year on the station
    fn = sprintf("%s/%sh%d.txt", station_id, station_id, i);
    if isfile(fn) % Check if the data file exists or not
        continue
    end

    try
        url = sprintf("https://www.ndbc.noaa.gov/view_text_file.php?filename=%sh%d.txt.gz&dir=data/historical/stdmet/", station_id, i);
        request = webread(url);
    catch
        fprintf("Failed to download buoy data for year %d. Most likely, the network is not connected or the file does not exist.", i)
        continue
    end

    fid = fopen(fn, 'w');           % Open the data file for writing data
    fprintf(fid, "%s", request);    % Write data to the file
    fclose(fid);                    % Close the data file
end

% Download current buoy data - CREATES BROKEN DATA FILE
% try
%     url = sprintf("https://www.ndbc.noaa.gov/data/realtime2/%s.txt", station_id);
%     request = webread(url);
% catch
%     disp("Failed to download buoy data. Check data connection!")
% end

% fn = sprintf("%s/%sh%d.txt", station_id, station_id, end_year+1);
% fid = fopen(fn, 'w');
% fprintf(fid, "%s", request);
% fclose(fid);

%% Import Data
fn = sprintf("%s/%s_data_%d_%d.csv", station_id, station_id, 2005, end_year);
if ~isfile(fn) % If the data table is not saved, then process all of the data and save it
    buoy_data = table();
    for i=start_year:1998 % For every year between the start year and 1998 - In 1999, the file format changed
        disp(i)
        buoy_data = vertcat(buoy_data, import_old_old_ndbc_data(sprintf("%s/%sh%02d.txt",station_id, station_id, i))); % Combine the buoy table together
    end
    for i=1999:2004 % For every year between 1999 and 2004 - In 2005, the file format changed
        disp(i)
        buoy_data = vertcat(buoy_data, import_old_ndbc_data(sprintf("%s/%sh%02d.txt",station_id, station_id, i))); % Combine the buoy table together
    end
    for i=2005:end_year % For every year between 2005 and end year
        disp(i)
        buoy_data = vertcat(buoy_data, import_ndbc_data(sprintf("%s/%sh%02d.txt", station_id, station_id, i))); % Combine the buoy tables together
    end
    buoy_data(ismember(buoy_data.WVHT,99.0),:)=[];
    buoy_data(ismember(buoy_data.WTMP,999.0),:)=[];
%     buoy_data = standardizeMissing(buoy_data, {999 99}); % Filter out unacceptable values
end

%% Display Surface and Histogram
% Plot wave surface and histogram and calculate statistics
figure
    subplot(3,1,1)
%     plot(buoy_data.DATE, buoy_data.WVHT) % Plot water levels over time
    plot(buoy_data.Time, buoy_data.WVHT) % Plot water levels over time
    title("Raw Wave Heights")
    ylabel("Wave Heights [m]")
    subplot(3,1,2)
    nbins = 50;
    histogram(buoy_data.WVHT, nbins); % Plot
    title("Wave Height Histogram (Bins = " + num2str(nbins) + ")")
    xlabel("Wave Heights [m]")
    ylabel("# of Occurances")
    subplot(3,1,3)
    plot(buoy_data.Time, buoy_data.WTMP) % Plot water temperature over time
    title("Water Temperature Over Time")
    ylabel("Water Temperature [°C]")
    
    saveas(gcf, sprintf("%s/wave_heights_temps_%s.png", station_id, station_id))
H_avg = mean(buoy_data.WVHT);
H_stdev = std(buoy_data.WVHT);
N = length(buoy_data.WVHT);
temp = sort(buoy_data.WVHT);
H_third = mean(temp(round(2/3*length(buoy_data.WVHT)):end));

%% PDF and CDF
% Plot the wave height probability distribution function (PDF) and cumulative 
% distribution function (CDF) for the provided data set.
figure
    hold on
    [count, edges] = histcounts(buoy_data.WVHT, nbins);
    for i=1:length(edges)-1 % For every bin edge value
        x(i) = (edges(i+1) + edges(i)) / 2; % Determine the bin center from the average of the edges
    end
    PDF = count / N;
    CDF = cumsum(PDF);
    [count, edges] = histcounts(buoy_data.WVHT, nbins, 'Normalization', 'pdf');
    plot(x, 100*PDF)
    plot(x, 100*CDF)
    plot(x, 100*count)
    title("Wave Height Probability Distributions")
    xlabel("Wave Height [m]")
    ylabel("Probability")
    ytickformat("percentage")
    legend("PDF (From CDF)", "CDF", "PDF (Historical)")
    hold off
    
    saveas(gcf, sprintf("%s/pdf_cdf_%s.png", station_id, station_id))

%% Estimate RI Wave Conditions
% Using the graphical historic method estimate the 50, 100 and 500 year 
% design wave conditions (Plot ln(Hn) vs. ln(n) where n = 1/(1-F(Hn)).

n = 1./(1-CDF);
n50 = (N/N_years)*50;
n100 = (N/N_years)*100;
n500 = (N/N_years)*500;

figure
    hold on
    plot(log(x),log(n))
    plot(x, log(n))
    yline(log(n50))
    yline(log(n100))
    yline(log(n500))
    title("Graphical Historical Estimate")
    xlabel("Wave Height (ln(H), H) [m]")
    ylabel("Return Interval Occurance (ln(n))")
    legend("Historical Data (ln)", "Historical Data", "50-year", "100-year", "500-year", 'Location', 'southeast')
    xlim([0 7])
    hold off
    
    saveas(gcf, sprintf("%s/graphical_historical_estimate_%s.png", station_id, station_id))

%% Extreme Analysis Distributions
% Which of the standard extreme analysis distributions (such as those 
% provided in the CEM Part II Chapter 8, Figure II-8-1 and 2) does the data
% best match: Rayleigh, FT-I (Gumbel), or Weibull (k=1, 1.4, and 2)?  
% (Plot all PDF’s and CDF’s including those from part a. together on one 
% figure and use the figure to arrive at your conclusion)

H_avg = mean(buoy_data.WVHT);
x1 = 0:0.05:14;

% Rayleigh Distribution Calculations
alpha1 = (0.886/H_avg)^2;
alpha2 = (0.463/std(buoy_data.WVHT))^2;
RayAlpha = mean([alpha1, alpha2]);
RayPDF = 2.*RayAlpha.*x1.*exp(-RayAlpha.*x1.^2);
RayCDF = 1-exp(-RayAlpha.*x1.^2);

% FT-I (Gumbell) Calculations
A = 0.779*std(buoy_data.WVHT);
B = H_avg - 0.45*std(buoy_data.WVHT);
GumCDF = exp(-exp(-(x1-B)/A));
GumPDF = diff(GumCDF) ./ diff(x1);

% Weibull (k=1.0)
A = std(buoy_data.WVHT);
B = H_avg - std(buoy_data.WVHT);
WeiCDF_1 = 1-exp(-((x1-B)/A).^1);
WeiPDF_1 = diff(WeiCDF_1) ./ diff(x1);

% Weibull (k=1.4)
A = 1.515*std(buoy_data.WVHT);
B = H_avg - 1.38*std(buoy_data.WVHT);
WeiCDF_1P4 = 1-exp(-((x1-B)/A).^1.4);
WeiPDF_1P4 = diff(WeiCDF_1P4) ./ diff(x1);

% Weibull (k=2.0)
A = 2.160*std(buoy_data.WVHT);
B = H_avg - 1.914*std(buoy_data.WVHT);
WeiCDF_2 = 1-exp(-((x1-B)/A).^2);
WeiPDF_2 = diff(WeiCDF_2) ./ diff(x1);

% Plotting
figure
    hold on
    plot(x, count*100, x, CDF*100)
    plot(x1, RayPDF*100, x1, RayCDF*100)
    plot(x1(2:end), GumPDF*100, x1, GumCDF*100)
    plot(x1(2:end), WeiPDF_1*100, x1, WeiCDF_1*100)
    plot(x1(2:end), WeiPDF_1P4*100, x1, WeiCDF_1P4*100)
    plot(x1(2:end), WeiPDF_2*100, x1, WeiCDF_2*100)
    ylim([0 100])
    title("Extreme Analysis Distributions")
    xlabel("Wave Height [m]")
    ylabel("Probability")
    ytickformat("percentage")
    legend("Historic PDF", "Historic CDF", ...
            "Rayleigh PDF", "Rayleigh PDF", ...
            "Gumbell PDF", "Gumbell CDF", ...
            "Weibull PDF (k=1.0)", "Weibull CDF (k=1.0)", ...
            "Weibull PDF (k=1.4)", "Weibull CDF (k=1.4)", ...
            "Weibull PDF (k=2.0)", "Weibull CDF (k=2.0)", ...
            "Location", "east")
    hold off
    saveas(gcf, sprintf("%s/extreme_distribution_%s.png", station_id, station_id))
    
%% Design Wave Estimation
% Using the Distributions in part c., estimate the 100 year and 500-year 
% return interval design wave. Plot a figure of return interval on the 
% x-axis and wave height on the y axis.

t = N_years/N; % Time of each sample in years.

% Using equation II-8-1 from CEM to find return intervals for each
Tr = t./(1-GumCDF);
Tr2 = t./(1-WeiCDF_1);
Tr3 = t./(1-WeiCDF_1P4);
Tr4 = t./(1-WeiCDF_2);
Tr5 = t./(1-RayCDF);
Tr6 = t./(1-CDF);

figure
    hold on
    plot(Tr, x1, 'm') 
    plot(Tr2, x1, 'c')
    plot(Tr3, x1, 'r')
    plot(Tr4, x1, 'g')
    plot(Tr5, x1, 'k')
    plot(Tr6, x(1:length(CDF)), 'b')
    set(gca, 'XScale', 'log')
    ylim([0 16])
    xlim([1 700])
    xline(50,'--k')
    xline(100,'--b')
    xline(500,'--r')
    grid on
    legend('Gumbell CDF', 'Weibull k=1', 'Weibull k=1.4', 'Weibull k=2', 'Rayleigh CDF', 'Historical', 'Location','northwest')
    title('Probability Distribution for Historic Wave Data')
    xlabel ('Return Interval in Years')
    ylabel ('Wave Height (m)')
    hold off
    saveas(gcf, sprintf("%s/probability_distribution_%s.png", station_id, station_id))

% Display
disp("The 50-year return interval wave is 9m (Weibull [k=1])")
disp("The 100-year return interval wave is 9.5m (Weibull [k=1])")