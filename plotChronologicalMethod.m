function plotChronologicalMethod( timeYearly, electricityNormalizedYearly, sortVectorMonthly, ...
                                                           sortVector, month, dStart, dEnd, CF, yearFlag)
%plotChronologicalMethodMonthly plots a nice graph showing which
%electricity use points are used for chronological method capacity factor
%calculation
figure(998); hold on;
CFpoint = round(CF*(dEnd-dStart));
title('Monthly chronological method')
plot(timeYearly(dStart:dEnd)/24,electricityNormalizedYearly(dStart:dEnd),'k'); hold on;
plot(timeYearly(dStart:dEnd)/24, electricityNormalizedYearly(sortVectorMonthly(month,1:(dEnd-dStart+1))),'.r','lineWidth',2)
plot(timeYearly(sortVectorMonthly(month,1:CFpoint))/24, ...
    electricityNormalizedYearly(sortVectorMonthly(month,1:CFpoint)),'^r');
datetick('x','mmm')
ylabel('normalized power');

if yearFlag
    figure(997); hold on;
    CFpoint = round(CF*length(timeYearly));
    title('Yearly chronological method')
    plot(timeYearly/24,electricityNormalizedYearly,'k'); hold on;
    plot(timeYearly/24, electricityNormalizedYearly(sortVector),'.r','lineWidth',2)
    plot(timeYearly(sortVector(1:CFpoint))/24, ...
        electricityNormalizedYearly(sortVector(1:CFpoint)),'^r');
    datetick('x','mmm')
    ylabel('normalized power');
    axis([0 timeYearly(end)/24 0.3 1])
    legend('normalized electricity consumption','sorted electricity consumption','top 30% of electricity consumption','Location','South')
    
    figure(999);
    subplot(121);
    startSummer = round(226.1*24*6);
    endSummer = round(227.2*24*6);
    CFpoint = round(CF*(endSummer-startSummer));
    dailySummerElectricity = electricityNormalizedYearly(startSummer:endSummer);
    [~,dailySortVector]=sort(dailySummerElectricity,'descend');
    plot(timeYearly(startSummer:endSummer)/24,dailySummerElectricity,'k'); hold on;
    plot(timeYearly(startSummer:endSummer)/24,dailySummerElectricity(dailySortVector),'.r','lineWidth',2)
    plot(timeYearly(startSummer + dailySortVector(1:CFpoint))/24, ...
        dailySummerElectricity(dailySortVector(1:CFpoint)),'^r');
    axis tight
    datetick('x','HH:MM','keeplimits')
    ylabel('normalized power');
    title('typical summer day')
    subplot(122);
    startWinter = round(25.1*24*6);
    endWinter = round(26.2*24*6);
    CFpoint = round(CF*(endWinter-startWinter));
    dailyWinterElectricity = electricityNormalizedYearly(startWinter:endWinter);
    [~,dailySortVector]=sort(dailyWinterElectricity,'descend');
    plot(timeYearly(startWinter:endWinter)/24,dailyWinterElectricity,'k'); hold on;
    plot(timeYearly(startWinter:endWinter)/24,dailyWinterElectricity(dailySortVector),'.r','lineWidth',2)
    plot(timeYearly(startWinter + dailySortVector(1:CFpoint))/24, ...
        dailyWinterElectricity(dailySortVector(1:CFpoint)),'^r');
    axis tight
    datetick('x','HH:MM','keeplimits')
    ylabel('normalized power');
    title('typical Winter day')
end
end

