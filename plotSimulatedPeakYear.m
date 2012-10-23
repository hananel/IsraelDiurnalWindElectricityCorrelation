function plotSimulatedPeakYear( maxDailyElectricity, dailyPeakCapacityFactor,yearlyCapacityFactor, stationName, year, stations, weights)
%plotSimulatedPeakYear 
global joker
figure(10010+year+joker); clf; hold on;
electricityInstalledPower = 15000; %[Mw]
maxDailyElectricity = maxDailyElectricity/max(maxDailyElectricity)*electricityInstalledPower; % renormalizing to electricityInstalledPower
% find installed power inorder to produce at least 400 Mw according to
% some rule
windInstalledPower = 1500; %400 / nanmean(dailyPeakCapacityFactor); %[Mw]

maxDailyWind = windInstalledPower*dailyPeakCapacityFactor;
maxDailyElectricityWithWind = maxDailyElectricity - maxDailyWind;

[ax,h1,h2] = plotyy(0.5:364.5,maxDailyElectricity,0.5:364.5, maxDailyWind);
hold on;
h3 = plot(0.5:364.5, maxDailyElectricityWithWind,'m','lineWidth',2);
datetick(ax(1), 'x','mmm'); datetick(ax(2), 'x','mmm'); 
set(ax(2),'Ylim',[0 3500])
set(ax(2),'YTick',[0 500 1000 1500])
set(ax(2),'YTickLabel',[0 500 1000 1500])

set(ax(1),'Ylim',[2500 15000])
set(ax(1),'YTick',[7500 10000 12500 15000])
set(ax(1),'YTickLabel', [7500 10000 12500 15000])

set(h1,'LineWidth',2)
set(h2,'LineWidth',2)
set(h1,'Color','k')
set(h2,'Color','b')

set(ax(2),'YColor','b')
set(ax(1),'YColor','k')

grid on;
grid(ax(2),'on')

% show Mw in normal numbers, not scientific
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%d'))

title({['simulation with ' stationName, ' data for ',num2str(year)],['Total installed capacity = ' num2str(electricityInstalledPower) ' Mw'], ...
       ['Total wind capacity = ', num2str(ceil(windInstalledPower)),' Mw']});
ylabel('Mw - electricity demand'); ylabel(ax(2),'Mw - wind')
c = {'Electricity cosumption','Wind power production','Reduced consumption with added wind capacity'};
h = [h1,h2,h3];
legend(h,c,'Location','southOutside')

maxNormal = electricityInstalledPower;
maxWind = max(maxDailyElectricityWithWind);
addedBenefit = round(maxNormal - maxWind);
plot([0.5 364.5],[maxNormal maxNormal],'k')
plot([0.5 364.5],[maxWind maxWind],'b--','linewidth',2)
text(20,maxWind*0.99,{[num2str(addedBenefit) ' Mw difference'], ...
                      ['Summer peak capacity credit = ',num2str(round(100*addedBenefit/windInstalledPower)),'%'], ...
                      ['Yearly capacity factor = ', num2str(round(100*yearlyCapacityFactor)), '%']})
end

