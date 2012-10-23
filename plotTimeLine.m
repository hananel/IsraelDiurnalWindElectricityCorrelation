function plotTimeLine(stationName,num,year,yearLegend,CF,timeYearly,electricityNormalizedYearly, ...
    stationWindEnergyNormalizedYearly,dailyPeakCapacityFactor,peakDailyTime,Prated,resultsDirectory)
%plotTimeLine plots a wind measuring stations energy production time line
%for a specific year vs. the electricity demand

figure(200+num*7+year);
set(gcf,'Units','normalized')
set(gcf, 'Position', [0.0304    0.4889    0.3615    0.3833])
set(gcf,'name',[stationName, ' time line for ', num2str(yearLegend)])
noNanPoints = ~isnan(stationWindEnergyNormalizedYearly);
timeVec = timeYearly(noNanPoints)/24;
windVec = stationWindEnergyNormalizedYearly(noNanPoints);
a = area(timeVec,windVec); hold on;
set(a,'FaceColor','g')
set(a,'EdgeColor','g')
plot(timeYearly/24,electricityNormalizedYearly,'k');
% because of the time it takes to plot the time line, showing a specific day makes sense
[~,showDay] = max(electricityNormalizedYearly(noNanPoints)); 
showDay = floor((showDay+find(noNanPoints>0,1))/24/6);
% finding maximum consumption in this day
loc = and(timeYearly/24>showDay,timeYearly/24<showDay+1);
[~,maxDailyElectricityLoc] = max(electricityNormalizedYearly(loc));
a = area([(timeYearly(maxDailyElectricityLoc+1)-2)/24+showDay, (timeYearly(maxDailyElectricityLoc+1)+2)/24+showDay],[Prated Prated]); hold on;  
set(a,'FaceColor','c')
set(a,'EdgeColor','c')
set(get(a,'Children'),'FaceAlpha',0.5)
title([stationName ' wind timeline for ' num2str(yearLegend) ' on maximum consumption day' ])
axis([showDay showDay+1 0 1])
datetick('x','dd.mmm HH:MM','keeplimits'); grid on;
zoomAdaptiveDateTicks('on')
print([resultsDirectory 'DiurnalPeakDay_', num2str(yearLegend), '_', strrep(stationName,' ',''),'.png'],'-dpng');
saveas(gcf, [resultsDirectory, 'DiurnalPeakDay_', num2str(yearLegend), '_', strrep(stationName,' ','')], 'fig')

figure(1600+num*7+year);
hold on; plot(timeVec,windVec,'r')
[h_ax,h_l] = plot2ylinscales(peakDailyTime,dailyPeakCapacityFactor,1/CF);
ylabel(h_ax(1),'CC')
ylabel(h_ax(2),'CC/CF')
set(h_ax(2),'Xtick',[])
set(h_l,'Marker','*'); set(h_l,'Color','k'); set(h_l,'LineStyle','none')
title({[stationName ' wind timeline for ' num2str(yearLegend) ' with peak capacity credit at each day'],['Yearly capacity factor = ' num2str(100*CF,2) '%']})
legend('wind','wind at daily peak','Location','Best')
axis tight
datetick('x','mmm','keeplimits'); grid on;
zoomAdaptiveDateTicks('on')
print([resultsDirectory 'dailyPeaks_', num2str(yearLegend), '_',strrep(stationName,' ',''),'.png'],'-dpng');
saveas(gcf, [resultsDirectory, 'dailyPeaks_', num2str(yearLegend), '_',  strrep(stationName,' ','')], 'fig')
end

