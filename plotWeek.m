function plotWeek( timeYearly, stationWindEnergyNormalizedYearlyInterp, electricityNormalizedYearly...
                                    , stationName, year, fig1, fig2, NIPnormalizedYearlyInterp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% typical summer week
figure(fig1); clf; hold on;
loc = and(timeYearly/24>30*7+7,timeYearly/24<30*7+14);
plot(timeYearly(loc)/24,stationWindEnergyNormalizedYearlyInterp(loc),'b')
title({'simulated wind energy - for a typical summer week',['station data used: ' stationName],num2str(year)})
datetick('x','dd mmm','keeplimits')
grid on
ylabel('normalized power')
if nargin>7
   plot(timeYearly(loc)/24,NIPnormalizedYearlyInterp(loc),'r') 
end
plot(timeYearly(loc)/24,electricityNormalizedYearly(loc),'k')
hold off;
if nargin>7
    legend('wind','solar','electricity')
else
    legend('wind','electricity')
end
figure(fig2); clf; hold on;
loc = and(timeYearly/24>30*1,timeYearly/24<30*1+7);
plot(timeYearly(loc)/24,stationWindEnergyNormalizedYearlyInterp(loc),'b')
title({'simulated wind energy - for a typical winter week',['station data used: ' stationName],num2str(year)})
datetick('x','dd mmm','keeplimits')
grid on
ylabel('normalized power')
if nargin>7
   plot(timeYearly(loc)/24,NIPnormalizedYearlyInterp(loc),'r') 
end
plot(timeYearly(loc)/24,electricityNormalizedYearly(loc),'k')
if nargin>7
    legend('wind','solar','electricity')
else
    legend('wind','electricity')
end
hold off;
end

