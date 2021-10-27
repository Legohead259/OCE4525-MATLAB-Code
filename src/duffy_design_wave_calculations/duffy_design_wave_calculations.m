% OCE4525 - Design Wave Calculations
% Braidan Duffy in collaboration with Gabriella Santiago and Clay Kenney
% October 8, 2021

close all; clear all;
load('TakeHome_dataset-1.mat')

% (100 points) Given the dataset provided (Historic 6 hour wave data), 
% determine the 50, 100 and 500 year design wave conditions.  
% Show all work, use MATLAB.  Provide your answers in this MS WORD Document.

%    A. (25pts) Plot the wave height probability distribution function (PDF) 
%               and cumulative distribution function (CDF) for the provided data set.
%    B. (20pts) Using the graphical historic method estimate the 50, 100 
%               and 500 year design wave conditions (Plot ln(Hn) vs. ln(n) where n = 1/(1-F(Hn)).
%    C. (30 pts) Which of the standard extreme analysis distributions 
%               (such as those provided in the CEM Part II Chapter 8, 
%               Figure II-8-1 and 2) does the data best match: Rayleigh, 
%               FT-I (Gumbel), or Weibull (k=1, 1.4, and 2)?  
%               (Plot all PDF’s and CDF’s including those from part a. 
%               together on one figure and use the figure to arrive at your conclusion)
%    D. (25 pts) Using the Distributions in part c., 
%               estimate the 100 year and 500 year return interval design wave. 
%               Plot a figure of return interval on the x axis and wave height on the y axis.  
%               (Hint: Derive equation using Eq. II-8-1 and II-8-2 along 
%               with the CDF equation from the table.  Solve II-8-2 for H in terms of Tr)

%% Part 0
% Plot wave surface and histogram and calculate statistics
figure
    subplot(2,1,1)
    plot(H) % Plot water levels over time
    title("Raw Wave Heights")
    ylabel("Wave Heights [m]")
    subplot(2,1,2)
    nbins = 50;
    histogram(H, nbins); % Plot
    title("Wave Height Histogram (Bins = " + num2str(nbins) + ")")
    xlabel("Wave Heights [m]")
    ylabel("# of Occurances")
   
H_avg = mean(H);
H_stdev = std(H);

%% Part A:
% Plot the wave height probability distribution function (PDF) and cumulative 
% distribution function (CDF) for the provided data set.
figure
    hold on
    [count, edges] = histcounts(H, nbins);
    for i=1:length(edges)-1 % For every bin edge value
        x(i) = (edges(i+1) + edges(i)) / 2; % Determine the bin center from the average of the edges
    end
    PDF = count / N;
    CDF = cumsum(PDF);
    [count, edges] = histcounts(H, nbins, 'Normalization', 'pdf');
    plot(x, 100*PDF)
    plot(x, 100*CDF)
    plot(x, 100*count)
    title("Wave Height Probability Distributions")
    xlabel("Wave Height [m]")
    ylabel("Probability")
    ytickformat("percentage")
    legend("PDF (From CDF)", "CDF", "PDF (Historical)")
    hold off

%% Part B:
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

%% Part C
% Which of the standard extreme analysis distributions (such as those 
% provided in the CEM Part II Chapter 8, Figure II-8-1 and 2) does the data
% best match: Rayleigh, FT-I (Gumbel), or Weibull (k=1, 1.4, and 2)?  
% (Plot all PDF’s and CDF’s including those from part a. together on one 
% figure and use the figure to arrive at your conclusion)

H_avg = mean(H);
x1 = 0:0.05:14;

% Rayleigh Distribution Calculations
alpha1 = (0.886/H_avg)^2;
alpha2 = (0.463/std(H))^2;
RayAlpha = mean([alpha1, alpha2]);
RayPDF = 2.*RayAlpha.*x1.*exp(-RayAlpha.*x1.^2);
RayCDF = 1-exp(-RayAlpha.*x1.^2);

% FT-I (Gumbell) Calculations
A = 0.779*std(H);
B = H_avg - 0.45*std(H);
GumCDF = exp(-exp(-(x1-B)/A));
GumPDF = diff(GumCDF) ./ diff(x1);

% Weibull (k=1.0)
A = std(H);
B = H_avg - std(H);
WeiCDF_1 = 1-exp(-((x1-B)/A).^1);
WeiPDF_1 = diff(WeiCDF_1) ./ diff(x1);

% Weibull (k=1.4)
A = 1.515*std(H);
B = H_avg - 1.38*std(H);
WeiCDF_1P4 = 1-exp(-((x1-B)/A).^1.4);
WeiPDF_1P4 = diff(WeiCDF_1P4) ./ diff(x1);

% Weibull (k=2.0)
A = 2.160*std(H);
B = H_avg - 1.914*std(H);
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
    
%% Part D
% Using the Distributions in part c., estimate the 100 year and 500-year 
% return interval design wave. Plot a figure of return interval on the 
% x-axis and wave height on the y axis. (Hint: Derive equation using 
% Eq. II-8-1 and II-8-2 along with the CDF equation from the table.  
% Solve II-8-2 for H in terms of Tr)

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

    % Display
    disp("The 50-year return interval wave is 7m (Rayleigh)")
    disp("The 100-year return interval wave is 7.25m (Rayleigh)")
    disp("The 500-year return interval wave is 7.7m (Rayleigh)")