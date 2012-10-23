function plotDailyWindSolarElectricityVisualization(electricityNormalizedDaily,stationWindEnergyNormalizedDaily,NIPDaily, monthsNum)
% produce a daily profile - which is an interannual average of a month -
% for solar, wind and electricity. for the front of the presentation :)

highNum  =1234;
months = {'January','February','March','April','May','June','July','August','September','October','November','December'};

% plotting for each month, or for the selected months only
if nargin<4
    monthsNum = 1:12;
end

% averaging all of the wind energy daily data
windEnergyDailyAveraged = nanmean(stationWindEnergyNormalizedDaily,3);

for month=monthsNum
    % normalizing solar and wind data
    NIPDaily(month,:) = NIPDaily(month,:)/max(NIPDaily(month,:));
    windEnergyDailyAveraged(month,:) = windEnergyDailyAveraged(month,:)/max(windEnergyDailyAveraged(month,:));
    % summing and renormalizing wind+solar
    renewableEnergyDaily = NIPDaily+windEnergyDailyAveraged;
    renewableEnergyDaily(month,:) = renewableEnergyDaily(month,:)/max(renewableEnergyDaily(month,:));
    % plotting
    figure(highNum+month); hold on;
    plot((0.5:23.5)/24,electricityNormalizedDaily(month,:),'k','lineWidth',3)
    plot((0.5:144)/6/24,windEnergyDailyAveraged(month,:),'b','lineWidth',2)
    plot((0.5:144)/6/24,NIPDaily(month,:),'r','lineWidth',2)
    plot((0.5:144)/6/24,renewableEnergyDaily(month,:),'g','lineWidth',2)
    legend('electricity consumption','wind power','solar power','wind+solar','Location','Best')
    title(['Normalized power trend for Israel - ' months{month}])
    xlabel('hour');
    ylabel('normalized power')
    axis tight
    datetick('x','HH:MM');
    set(gca,'YTick',[0 1])
    grid on;
end