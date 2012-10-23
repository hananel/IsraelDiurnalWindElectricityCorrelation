function plotIEC( tYearly,electricityNormalizedYearly, electricityNormalizedYearlySorted, resultsDirectory)
%plotIEC plots sorted and unsorted IEC data
monthString = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
figure(15);
subplot(1,2,1); set(gca,'fontsize',18)
plot(tYearly/24,electricityNormalizedYearly,'k','LineWidth',3);
datetick('x','mmm','keeplimits')
set(gca,'ytick',[0 1])
xlabel('month'); ylabel('normalized load');
title('typical yearly load');
subplot(1,2,2)
plot((1:length(electricityNormalizedYearlySorted))/24/6/365*100,electricityNormalizedYearlySorted,'k','LineWidth',3);
set(gca,'xtick',[0 30 50 100]);
set(gca,'ytick',[0 1]);
grid on;
axis([0 100 0 1])
xlabel('percent of maximum load'); ylabel('normalized load');
title('sorted yearly load');
print([resultsDirectory 'typicalLoadIEC_2011_2030.png'],'-dpng');  

end

