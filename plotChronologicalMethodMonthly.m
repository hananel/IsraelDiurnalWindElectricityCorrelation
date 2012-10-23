function plotChronologicalMethod( timeYearly, electricityNormalizedYearly, sortVectorMonthly, ...
                                                           sortVector, month, dStart, dEnd, CF)
%plotChronologicalMethodMonthly plots a nice graph showing which
%electricity use points are used for chronological method capacity factor
%calculation
figure(997); hold on;
CFpoint = round(CF*(dEnd-dStart));
title('Monthly chronological method')
plot(timeYearly(dStart:dEnd)/24,electricityNormalizedYearly(dStart:dEnd),'k'); hold on;
plot(timeYearly(dStart:dEnd)/24, electricityNormalizedYearly(sortVectorMonthly(month,1:(dEnd-dStart+1))),'.r','lineWidth',2)
plot(timeYearly(sortVectorMonthly(month,1:CFpoint))/24, ...
    electricityNormalizedYearly(sortVectorMonthly(month,1:CFpoint)),'^r');
datetick('x','mmm')
ylabel('normalized power');

figure(998); hold on;
title('Yearly chronological method')
plot(timeYearly/24,electricityNormalizedYearly,'k'); hold on;
plot(timeYearly/24, electricityNormalizedYearly(sortVector),'.r','lineWidth',2)
plot(timeYearly(sortVector)/24, ...
    electricityNormalizedYearly(sortVector),'^r');
datetick('x','mmm')
ylabel('normalized power');

figure(999);
title('Daily peak capacity method')

end

