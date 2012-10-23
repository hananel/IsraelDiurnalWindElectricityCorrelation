function plotChronologicalMethodMonthly( timeYearly, electricityNormalizedYearly, electricityNormalizedYearly,sortVectorMonthly, ...
                                                           month, dStart, dEnd,  )
%plotChronologicalMethodMonthly plots a nice graph showing which
%electricity use points are used for chronological method capacity factor
%calculation
figure(999); hold on;
CFpoint = round(CF*(dEnd-dStart));
plot(timeYearly(dStart:dEnd),electricityNormalizedYearly(dStart:dEnd),'k'); hold on;
plot(timeYearly(dStart:dEnd), electricityNormalizedYearly(sortVectorMonthly(month,dStart:dEnd)),'.r','lineWidth',4)
plot(timeYearly(sortVectorMonthly(month,dStart:CFpoint)), ...
    electricityNormalizedYearly(sortVectorMonthly(month,dStart:CFpoint)),'^g');
plot([CFpoint CFpoint],[min(electricityNormalizedYearly) max(electricityNormalizedYearly)],'r','lineWidth',2)
end

