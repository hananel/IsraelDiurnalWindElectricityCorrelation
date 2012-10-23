function plotIsraelYear(year, timeYearly, EWnormYIsraelSort, ...
    NIPnormYSort, EWNIPYSort, EWnormYIsrael, NIPnormY, EWNIPY, ...
    electricityNormalizedYearly)
%plotIsraelYear plots data from aggregate of IMS stations

figure(3000+year); hold on; 
set(gcf,'name','Israel sorted time line')
plot(timeYearly/24,EWnormYIsraelSort,'b'); hold on
plot(timeYearly/24,NIPnormYSort,'r')
plot(timeYearly/24,EWNIPYSort,'g')
legend('wind','solar','both','Location','SouthEast')
title(['Israel sorted timeline for ' num2str(year)])

figure(4000+year); hold on; 
set(gcf,'name','Israel time line') 
grid on;
plot(timeYearly/24,EWnormYIsrael,'b');
datetick('x','mmm')
plot(timeYearly/24,NIPnormY,'r')
plot(timeYearly/24,EWNIPY,'g')
plot((1:length(electricityNormalizedYearly))/24/6,electricityNormalizedYearly,'k');
legend('wind','solar','wind+solar','electricity','Location','SouthEast')
title(['Israel timeline for ' num2str(year)])

end

