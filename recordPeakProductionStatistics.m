function [monthlyPeakCapacityFactorMinimum, monthlyPeakCapacityFactorMaximum,monthlyPeakCapacityFactorAverage, monthlyPeakCapacityFactorStd, ...
    yearlyPeakCapacityFactorMinimum, yearlyPeakCapacityFactorMaximum, yearlyPeakCapacityFactorAverage, yearlyPeakCapacityFactorStd, ...
    peakDailyTime, dailyPeakCapacityFactor, maxDailyElectricity, yearCounter] = ...
    recordPeakProductionStatistics(timeYearly, electricityNormalizedYearly, stationWindEnergyNormalizedYearlyInterp, ...
    Prated, yearCounter, windDataYearVec, Num, year, minNan, ...
    monthlyPeakCapacityFactorMinimum, monthlyPeakCapacityFactorMaximum, ...
    monthlyPeakCapacityFactorAverage, monthlyPeakCapacityFactorStd, ...
    yearlyPeakCapacityFactorMinimum, yearlyPeakCapacityFactorMaximum, ...
    yearlyPeakCapacityFactorAverage, yearlyPeakCapacityFactorStd, ...
    stationName)
%recordPeakProductionStatistics recording maximum,minimum,average and std of capacity factor
% at peak consumption for each station
% TODO - add a daily CC % number - to check a precentage around the maximum and
% not just the max 10 minute load
% which is equivalent to 10/(24*60) =~ 0.7%
% width - [hr] - 7 hours width around peak, equivalent to 30%. 12 hours 50%,
%and the original 10 minute width is equivalent to 0.7%

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
global debug; global report
global reportName
for i=0:364 % walking over days
    % finding maximum consumption in this day
    loc = and(timeYearly/24>i,timeYearly/24<i+1);
    temp = electricityNormalizedYearly(loc);
    [maxDailyElectricity(i+1),maxDailyElectricityLoc] = max(temp);
    peakDailyTime(i+1) = maxDailyElectricityLoc/6/24 + i;
    % deciding on width - according to Yoel-Hanan meeting 29/8/12
    % 0.4 Gw / 15 Gw = 2.67% of normalized consumption
    % becuse IEC work according to 0.4 Gw production unit blocks, and 15 Gw
    % is the total installed capacity in Israel in a few years.
    % TODO (perhaps) change delta according to time of year? smaller at
    % winter??
    delta = 0.4/15;
    
    % checking how many intersection lines we have (one or double peak)
    Xpoly = polyxpoly([0 143],[maxDailyElectricity(i+1)-delta maxDailyElectricity(i+1)-delta],0:length(temp)-1,temp);
    if length(Xpoly)>2
        if debug
            figure(994); hold on; plot(temp)
            plot(Xpoly,ones(1,length(Xpoly))*(maxDailyElectricity(i+1)-delta),'r-*')
        end
        % right side
        widthLeft1 = floor(Xpoly(end-1));
        widthRight1 = ceil(Xpoly(end));
        startLoc1 = i*24*6 + widthLeft1;
        endLoc1 = i*24*6 + widthRight1;
        % left side
        widthLeft2 = floor(Xpoly(1));
        widthRight2 = ceil(Xpoly(2));
        startLoc2 = i*24*6 + widthLeft2;
        endLoc2 = i*24*6 + widthRight2;
        dailyPeakCapacityFactor(i+1) = nanmin([stationWindEnergyNormalizedYearlyInterp(startLoc1:endLoc1), ...
            stationWindEnergyNormalizedYearlyInterp(startLoc2:endLoc2)])/Prated;
        widthLeft = widthLeft1; widthRight = widthRight1; startLoc = startLoc2; endLoc = endLoc1;
        dailyWidth(i+1) = (widthRight2-widthLeft2+widthRight1-widthLeft1)/6;
    else
        widthLeft = find(temp>(maxDailyElectricity(i+1)-delta),1,'first');
        widthRight = find(temp>(maxDailyElectricity(i+1)-delta),1,'last');
        % recording peak+/-width/2 capacity factor
        startLoc = i*24*6+widthLeft;
        endLoc = i*24*6+widthRight;
        dailyPeakCapacityFactor(i+1) = nanmin(stationWindEnergyNormalizedYearlyInterp(startLoc:endLoc))/Prated;
        dailyWidth(i+1) = (widthRight-widthLeft)/6;
    end
    if debug
        figure(995); hold on; plot(i+(1:length(electricityNormalizedYearly(loc)))/length(electricityNormalizedYearly(loc)),electricityNormalizedYearly(loc),'k', ...
            [startLoc endLoc]/144,[maxDailyElectricity(i+1)-delta maxDailyElectricity(i+1)-delta],'r')
        text(peakDailyTime(i+1),maxDailyElectricity(i+1)+0.02,num2str((widthRight-widthLeft)/6,2))
        title(['day ' num2str(i) ', capacity factor ', num2str(dailyPeakCapacityFactor(i+1))])
    end
end

if debug
    figure(993); subplot(212);
    plot(dailyWidth)
    subplot(211);
    plot(dailyPeakCapacityFactor,'*-b');
    title(stationName)
end

% calculating maximum/minimum/average yearly values
yearlyPeakCapacityFactorMinimum(Num,year-windDataYearVec(1)+1) = nanmin(dailyPeakCapacityFactor);
yearlyPeakCapacityFactorMaximum(Num,year-windDataYearVec(1)+1) = nanmax(dailyPeakCapacityFactor);
yearlyPeakCapacityFactorAverage(Num,year-windDataYearVec(1)+1) = nanmean(dailyPeakCapacityFactor);
yearlyPeakCapacityFactorStd(Num,year-windDataYearVec(1)+1) = nanstd(dailyPeakCapacityFactor);

% calculating maximum/minimum/average monthly values

for i=1:12
    dStart = datenum(0,i,0,0,0,0)+1;
    dEnd = min(datenum(0,i+1,0,0,0,0),365);
    nanNum = sum(isnan(dailyPeakCapacityFactor(dStart:dEnd)));
    if nanNum<=minNan
        yearCounter(i) = yearCounter(i) + 1;
        monthlyPeakCapacityFactorMinimum(Num,i) = nanmin([monthlyPeakCapacityFactorMinimum(Num,i) ...
            , dailyPeakCapacityFactor(dStart:dEnd)]);
        monthlyPeakCapacityFactorMaximum(Num,i) = nanmax([monthlyPeakCapacityFactorMaximum(Num,i) ...
            , dailyPeakCapacityFactor(dStart:dEnd)]);
        monthlyPeakCapacityFactorAverage(Num,i) = (monthlyPeakCapacityFactorAverage(Num,i).*(yearCounter(i)-1) ...
            + nanmean(dailyPeakCapacityFactor(dStart:dEnd)))/yearCounter(i);
        monthlyPeakCapacityFactorStd(Num,i)     = (monthlyPeakCapacityFactorStd(Num,i).*(yearCounter(i)-1) ...
            + nanstd(dailyPeakCapacityFactor(dStart:dEnd)))/yearCounter(i);
    end
    if debug
        noNan = ~(isnan(dailyPeakCapacityFactor(dStart:dEnd)));
        temp = dailyPeakCapacityFactor(dStart:dEnd);
        [yy,xx] = hist(temp(noNan));
        figure(991);
        subplot(4,3,i); bar(xx,yy/nanmean(dailyPeakCapacityFactor));
        xlabel('CC/CF'); ylabel('% occurance'); title(months{i});
        
    end
end
%suptitle({'Peak shaving capability:','daily peak capacity credit normalized by yearly capacity factor'});
if debug
    hold off;
end
end

