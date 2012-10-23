 % copy-paste from CapacityCreditWindIsrael unfortonately...
 
 % loading solar data
[~,~,~, ~, NIP, ~,~,ySolar,~,~] = loadIMSradiation(75,5003,windDataDirectory,0);
IsraelSolarDataYearVec = unique(ySolar);
[y, ~, ~, ~, ~, ~] = datevec(tVec);
IsraelWindDataYearVec = unique(y);
IsraelWindYearlyPeakCapacityFactorMinimum = NaN*ones(1,length(IsraelWindDataYearVec));
IsraelWindYearlyPeakCapacityFactorMaximum = NaN*zeros(1,length(IsraelWindDataYearVec));
IsraelWindYearlyPeakCapacityFactorAverage = zeros(1,length(IsraelWindDataYearVec));
IsraelWindYearlyPeakCapacityFactorStd = zeros(1,length(IsraelWindDataYearVec));
IsraelSolarYearlyPeakCapacityFactorMinimum = NaN*ones(1,length(IsraelSolarDataYearVec));
IsraelSolarYearlyPeakCapacityFactorMaximum = NaN*zeros(1,length(IsraelSolarDataYearVec));
IsraelSolarYearlyPeakCapacityFactorAverage = zeros(1,length(IsraelSolarDataYearVec));
IsraelSolarYearlyPeakCapacityFactorStd = zeros(1,length(IsraelSolarDataYearVec));
IsraelCombinedYearlyPeakCapacityFactorMinimum = NaN*ones(1,length(IsraelWindDataYearVec));
IsraelCombinedYearlyPeakCapacityFactorMaximum = NaN*zeros(1,length(IsraelWindDataYearVec));
IsraelCombinedYearlyPeakCapacityFactorAverage = zeros(1,length(IsraelWindDataYearVec));
IsraelCombinedYearlyPeakCapacityFactorStd = zeros(1,length(IsraelWindDataYearVec));
yearCounter = zeros(12,1);
if length(IsraelWindDataYearVec)>6 IsraelWindDataYearVec=2006:2011; end
for year=IsraelWindDataYearVec
    
    % WIND DATA
    temp = WindEnergyNormalizedAveraged(y==year);
    WindEnergyNormalizedAveragedYearly(1:length(temp),year-IsraelWindDataYearVec(1)+1) = temp;
    % interpolating to timeYearly
    % this works in octave tOneYearWind = t(find(y==year))-t(find(y==year))(1); %[day]
    tOneYearWind = tVec(y==year);
    tOneYearWind = tOneYearWind - tOneYearWind(1);
    WindEnergyNormalizedAveragedYearlyInterp = interp1(tOneYearWind, WindEnergyNormalizedAveragedYearly(1:length(temp), year - IsraelWindDataYearVec(1) + 1), timeYearly / 24);
    CFIsrael(year-IsraelWindDataYearVec(1)+1) = nansum(WindEnergyNormalizedAveragedYearlyInterp)/length(WindEnergyNormalizedAveragedYearlyInterp)/Prated;
    % arranging the wind data in the sorted order of the load
    WindEnergyNormalizedAveragedYearlyInterpSort = WindEnergyNormalizedAveragedYearlyInterp(sortVector);
    
    % 8.
    % SOLAR DATA: doing the same for solar data
    % normalizing (same as above - all the 6 year meteorological data)
    NIPnormalized = NIP/max(NIP)*Prated; % making an equivalent solar and wind installation
    temp = NIPnormalized(y==year);
    NIPnormalizedYearly(1:length(temp),year-windDataYearVec(1)+1) = temp;
    % interpolating to timeYearly
    tOneYearWind = t(find(y==year));
    tOneYearWind = tOneYearWind - tOneYearWind(1);
    NIPnormalizedYearlyInterp = interp1(tOneYearWind,NIPnormalizedYearly(1:length(temp),year-windDataYearVec(1)+1),timeYearly/24);
    % arranging the solar data in the sorted order of the load
    NIPnormalizedYearlyInterpSort = NIPnormalizedYearlyInterp(sortVector);
    
    % 9.
    % WIND+SOLAR DATA
    % YOEL - check the normalization rule. is it appropriate?
    WindEnergyNIPnormalizedYearlyInterp = NIPnormalizedYearlyInterp + WindEnergyNormalizedAveragedYearlyInterp;
    % renormalizing the sum of the solar and wind stations
    WindEnergyNIPnormalizedYearlyInterp = WindEnergyNIPnormalizedYearlyInterp/max(WindEnergyNIPnormalizedYearlyInterp)*Prated;
    % arranging the combined data in the sorted order of the load
    WindEnergyNIPnormalizedYearlyInterpSort = WindEnergyNIPnormalizedYearlyInterp(sortVector);
    
    % calculating CC
    percentCounter = 0;
    for percent=percentVec
        percentCounter = percentCounter+1;
        loc = round(24*365*percent/100*6); % the 6 factor is for the 1 hour resolution for the load data vs. the 10 minute wind data
        % Israel wind
        CFperIsrael(year-IsraelWindDataYearVec(1)+1,percentCounter) = nansum(WindEnergyNormalizedAveragedYearlyInterpSort(1:loc))/loc/Prated;
        % Israel solar
        CFperIsraelS(year-IsraelWindDataYearVec(1)+1,percentCounter) = nansum(NIPnormalizedYearlyInterpSort(1:loc))/loc/Prated;
        % Israel wind+solar (same rated power for both)
        CFperIsraelSW(year-IsraelWindDataYearVec(1)+1,percentCounter) = nansum(WindEnergyNIPnormalizedYearlyInterpSort(1:loc))/loc/Prated;
    end
    CFIsrael(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsrael(year-IsraelWindDataYearVec(1)+1,end));
    CFIsrael1(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsrael(year-IsraelWindDataYearVec(1)+1,loc1));
    CFIsrael30(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsrael(year-IsraelWindDataYearVec(1)+1,loc30));
    CFIsrael50(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsrael(year-IsraelWindDataYearVec(1)+1,loc50));
    CFIsraelS1(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsraelS(year-IsraelWindDataYearVec(1)+1,loc1));
    CFIsraelS30(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsraelS(year-IsraelWindDataYearVec(1)+1,loc30));
    CFIsraelS50(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsraelS(year-IsraelWindDataYearVec(1)+1,loc50));
    CFIsraelSW1(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsraelSW(year-IsraelWindDataYearVec(1)+1,loc1));
    CFIsraelSW30(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsraelSW(year-IsraelWindDataYearVec(1)+1,loc30));
    CFIsraelSW50(year-IsraelWindDataYearVec(1)+1) = nanmean(CFperIsraelSW(year-IsraelWindDataYearVec(1)+1,loc50));
    
    % recording maximum,minimum,average and std of capacity factor
    % at peak consumption for wind, solar, and combined
    
    % Wind
    [IsraelWindMonthlyPeakCapacityFactorMinimum, IsraelWindMonthlyPeakCapacityFactorMaximum,IsraelWindMonthlyPeakCapacityFactorAverage, IsraelWindMonthlyPeakCapacityFactorStd, ...
        IsraelWindYearlyPeakCapacityFactorMinimum,  IsraelWindYearlyPeakCapacityFactorMaximum, IsraelWindYearlyPeakCapacityFactorAverage,  IsraelWindYearlyPeakCapacityFactorStd, ...
        peakDailyTime, IsraelWindDailyPeakCapacityFactor, maxDailyElectricity] = ...
        recordPeakProductionStatistics(timeYearly, electricityNormalizedYearly, WindEnergyNormalizedAveragedYearlyInterp, ...
        Prated, yearCounter, IsraelWindDataYearVec, 1, year, minNan, ...
        IsraelWindMonthlyPeakCapacityFactorMinimum, IsraelWindMonthlyPeakCapacityFactorMaximum, ...
        IsraelWindMonthlyPeakCapacityFactorAverage, IsraelWindMonthlyPeakCapacityFactorStd, ...
        IsraelWindYearlyPeakCapacityFactorMinimum, IsraelWindYearlyPeakCapacityFactorMaximum, ...
        IsraelWindYearlyPeakCapacityFactorAverage, IsraelWindYearlyPeakCapacityFactorStd, 'Israel wind');
    % simulation of 15000Mw/1500Mw scenario
    joker = 1;
    plotSimulatedPeakYear( maxDailyElectricity, IsraelWindDailyPeakCapacityFactor, CFIsrael(year-IsraelWindDataYearVec(1)+1), simulationName, year)
    if year==2011
        plotWeek( timeYearly, WindEnergyNormalizedAveragedYearlyInterp, electricityNormalizedYearly...
            , 'Israel combined', year, 11003, 11004, NIPnormalizedYearlyInterp)
    end
    % Solar
    [IsraelSolarMonthlyPeakCapacityFactorMinimum, IsraelSolarMonthlyPeakCapacityFactorMaximum,IsraelSolarMonthlyPeakCapacityFactorAverage, IsraelSolarMonthlyPeakCapacityFactorStd, ...
        IsraelSolarYearlyPeakCapacityFactorMinimum,  IsraelSolarYearlyPeakCapacityFactorMaximum, IsraelSolarYearlyPeakCapacityFactorAverage,  IsraelSolarYearlyPeakCapacityFactorStd, ...
        peakDailyTime, IsraelSolarDailyPeakCapacityFactor] = ...
        recordPeakProductionStatistics(timeYearly, electricityNormalizedYearly, NIPnormalizedYearlyInterp, ...
        Prated, yearCounter, IsraelSolarDataYearVec, 1, year, minNan, ...
        IsraelSolarMonthlyPeakCapacityFactorMinimum, IsraelSolarMonthlyPeakCapacityFactorMaximum, ...
        IsraelSolarMonthlyPeakCapacityFactorAverage, IsraelSolarMonthlyPeakCapacityFactorStd, ...
        IsraelSolarYearlyPeakCapacityFactorMinimum, IsraelSolarYearlyPeakCapacityFactorMaximum, ...
        IsraelSolarYearlyPeakCapacityFactorAverage, IsraelSolarYearlyPeakCapacityFactorStd, 'Israel solar');
    
    % Combined
    [IsraelCombinedMonthlyPeakCapacityFactorMinimum, IsraelCombinedMonthlyPeakCapacityFactorMaximum,IsraelCombinedMonthlyPeakCapacityFactorAverage, IsraelCombinedMonthlyPeakCapacityFactorStd, ...
        IsraelCombinedYearlyPeakCapacityFactorMinimum,  IsraelCombinedYearlyPeakCapacityFactorMaximum, IsraelCombinedYearlyPeakCapacityFactorAverage,  IsraelCombinedYearlyPeakCapacityFactorStd, ...
        peakDailyTime, IsraelCombinedDailyPeakCapacityFactor] = ...
        recordPeakProductionStatistics(timeYearly, electricityNormalizedYearly, WindEnergyNIPnormalizedYearlyInterp, ...
        Prated, yearCounter, IsraelWindDataYearVec, 1, year, minNan, ...
        IsraelCombinedMonthlyPeakCapacityFactorMinimum, IsraelCombinedMonthlyPeakCapacityFactorMaximum, ...
        IsraelCombinedMonthlyPeakCapacityFactorAverage, IsraelCombinedMonthlyPeakCapacityFactorStd, ...
        IsraelCombinedYearlyPeakCapacityFactorMinimum, IsraelCombinedYearlyPeakCapacityFactorMaximum, ...
        IsraelCombinedYearlyPeakCapacityFactorAverage, IsraelCombinedYearlyPeakCapacityFactorStd, 'Israel combined');
    
    % plotting combined stations - in year loop
    if debug
        plotIsraelYear( year, timeYearly, ...
            WindEnergyNormalizedAveragedYearlyInterpSort, NIPnormalizedYearlyInterpSort, WindEnergyNIPnormalizedYearlyInterpSort, ...
            WindEnergyNormalizedAveragedYearlyInterp, NIPnormalizedYearlyInterp, WindEnergyNIPnormalizedYearlyInterp, ...
            electricityNormalizedYearly)
        
    end
    sp1 = floor(length(IsraelWindDataYearVec)/2);
    sp2 = floor(length(IsraelWindDataYearVec)/sp1);
    if sp1*sp2<length(IsraelWindDataYearVec)
        sp2 = sp2+1;
    end
    % plotCC(percentVec, sp1,sp2,Num,year,year-IsraelWindDataYearVec(1)+1, ...
    %    CFperIsrael(year-IsraelWindDataYearVec(1)+1,:),0, CFperIsraelS(year-IsraelWindDataYearVec(1)+1,:), CFperIsraelSW(year-IsraelWindDataYearVec(1)+1,:));
end
% suptitle('Israel scenario - capacity Credit for precentage of peak demand')
% set(gcf,'Units','normalized')
% set(gcf,'Position',[0.126647 -0.294271 1.13031 1.10417])
% legend('wind','solar','wind+solar','Location','Best')
% saveas(gcf, [resultsDirectory 'IsraelCC_'], 'fig')

% plotting combined stations

figure(400)
loc = find(IsraelWindDataYearVec>0);
bar((ones(3,1)*IsraelWindDataYearVec)',[CFIsrael50(loc);CFIsraelS50(loc);CFIsraelSW50(loc)]'); hold on;
legend('Wind','Solar','Wind+Solar','Location','North');
minCF50 = floor(min(min([CFIsrael50(loc);CFIsraelS50(loc);CFIsraelSW50(loc)]))*10)/10-0.03; maxCF50 = ceil(max(max([CFIsrael50(loc);CFIsraelS50(loc);CFIsraelSW50(loc)]))*100)/100+0.01;
[h_ax,h_l] = plot2ylinscales([min(IsraelWindDataYearVec),max(IsraelWindDataYearVec)],[minCF50 maxCF50],1/nanmean(CFIsrael(loc)));
axis([min(IsraelWindDataYearVec)-0.5,max(IsraelWindDataYearVec)+0.5,minCF50 maxCF50]);
ylabel(h_ax(1),'50% CC')
ylabel(h_ax(2),'CC/CF')
lin_fun = inline([num2str(1/nanmean(CFIsrael)) '*x']);
set(h_ax(2),'ylim',lin_fun(get(h_ax(1),'ylim')));
set(h_l,'LineStyle','none'); set(h_ax(2),'Xtick',[])
h1 = title({'Israel 50% CC',['Scenario: Ulimit = ' num2str(Ulimit,2) ' m/s, power density = ' num2str(normPower) ' watt/m^2'],['interannual wind plants capacity factor = ' num2str(100*nanmean(CFIsrael),2) '%']});
axpos = get(gca,'pos');
extent = get(h1,'extent');
set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-2.5*extent(4)])
saveas(gcf, [resultsDirectory 'Israel50CC'], 'fig')

figure(401)
bar((ones(3,1)*IsraelWindDataYearVec)',[CFIsrael30(loc);CFIsraelS30(loc);CFIsraelSW30(loc)]'); hold on;
legend('Wind','Solar','Wind+Solar','Location','North');
minCF30 = floor(min(min([CFIsrael30(loc);CFIsraelS30(loc);CFIsraelSW30(loc)]))*10)/10-0.03; maxCF30 = ceil(max(max([CFIsrael30(loc);CFIsraelS30(loc);CFIsraelSW30(loc)]))*100)/100+0.01;
[h_ax,h_l] = plot2ylinscales([min(IsraelWindDataYearVec),max(IsraelWindDataYearVec)],[minCF30 maxCF30],1/nanmean(CFIsrael(loc)));
axis([min(IsraelWindDataYearVec)-0.5,max(IsraelWindDataYearVec)+0.5,minCF30 maxCF30]);
ylabel(h_ax(1),'30% CC')
ylabel(h_ax(2),'CC/CF')
lin_fun = inline([num2str(1/nanmean(CFIsrael)) '*x']);
set(h_ax(2),'ylim',lin_fun(get(h_ax(1),'ylim')));
set(h_l,'LineStyle','none'); set(h_ax(2),'Xtick',[])
h2 = title({'Israel 30% CC',['Scenario: Ulimit = ' num2str(Ulimit,2) ' m/s, power density = ' num2str(normPower) ' watt/m^2'],['interannual wind plants capacity factor = ' num2str(100*nanmean(CFIsrael),2) '%']});
axpos = get(gca,'pos');
extent = get(h2,'extent');
set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-2.5*extent(4)])
saveas(gcf, [resultsDirectory 'mIsrael30CC'], 'fig')

figure(402)
bar((ones(3,1)*IsraelWindDataYearVec)',[CFIsrael1(loc);CFIsraelS1(loc);CFIsraelSW1(loc)]'); hold on;
legend('Wind','Solar','Wind+Solar','Location','North');
minCF1 = floor(min(min([CFIsrael1(loc);CFIsraelS1(loc);CFIsraelSW1(loc)]))*10)/10-0.03; maxCF1 = ceil(max(max([CFIsrael1(loc);CFIsraelS1(loc);CFIsraelSW1(loc)]))*100)/100+0.01;
[h_ax,h_l] = plot2ylinscales([min(IsraelWindDataYearVec),max(IsraelWindDataYearVec)],[minCF1 maxCF1],1/nanmean(CFIsrael(loc)));
axis([min(IsraelWindDataYearVec)-0.5,max(IsraelWindDataYearVec)+0.5,minCF1 maxCF1]);
ylabel(h_ax(1),'1% CC')
ylabel(h_ax(2),'CC/CF')
lin_fun = inline([num2str(1/nanmean(CFIsrael)) '*x']);
set(h_ax(2),'ylim',lin_fun(get(h_ax(1),'ylim')));
set(h_l,'LineStyle','none'); set(h_ax(2),'Xtick',[])
h3 = title({'Israel 1% CC',['Scenario: Ulimit = ' num2str(Ulimit,2) ' m/s, power density = ' num2str(normPower) ' watt/m^2'],['interannual wind plants capacity factor = ' num2str(100*nanmean(CFIsrael),2) '%']});
axpos = get(gca,'pos');
extent = get(h3,'extent');
set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-1*extent(4)])
saveas(gcf, [resultsDirectory 'Israel1CC'], 'fig')

% plotting statistics of daily peak wind energy capacity
disp(' plotting peak statistics for Israel wind')
plotPeak('Israel wind', 1, IsraelWindYearlyPeakCapacityFactorMinimum, IsraelWindYearlyPeakCapacityFactorMaximum, ...
    IsraelWindYearlyPeakCapacityFactorAverage, IsraelWindYearlyPeakCapacityFactorStd, ...
    IsraelWindMonthlyPeakCapacityFactorMinimum, IsraelWindMonthlyPeakCapacityFactorMaximum, ...
    IsraelWindMonthlyPeakCapacityFactorAverage, IsraelWindMonthlyPeakCapacityFactorStd, ...
    resultsDirectory, nanmean(CFIsrael),IsraelWindDataYearVec);
% plotting monthly distribution of Capacity Credit and Factor for
% comparison
%plotMonthlyCC('Israel wind', CFMonthlyIsrael, CF50MonthlyIsrael, CF30MonthlyIsrael, CF1MonthlyIsrael, nanmean(CFIsrael))

% plotting statistics of daily peak solar energy capacity
disp(' plotting peak statistics for Israel solar')
plotPeak('Israel solar', 1, IsraelSolarYearlyPeakCapacityFactorMinimum, IsraelSolarYearlyPeakCapacityFactorMaximum, ...
    IsraelSolarYearlyPeakCapacityFactorAverage, IsraelSolarYearlyPeakCapacityFactorStd, ...
    IsraelSolarMonthlyPeakCapacityFactorMinimum, IsraelSolarMonthlyPeakCapacityFactorMaximum, ...
    IsraelSolarMonthlyPeakCapacityFactorAverage, IsraelSolarMonthlyPeakCapacityFactorStd, ...
    resultsDirectory, nanmean(CFIsrael),IsraelSolarDataYearVec);

% plotting statistics of daily peak combined energy capacity
disp(' plotting peak statistics for Israel combined solar + wind')
plotPeak('Israel combined', 1, IsraelCombinedYearlyPeakCapacityFactorMinimum, IsraelCombinedYearlyPeakCapacityFactorMaximum, ...
    IsraelCombinedYearlyPeakCapacityFactorAverage, IsraelCombinedYearlyPeakCapacityFactorStd, ...
    IsraelCombinedMonthlyPeakCapacityFactorMinimum, IsraelCombinedMonthlyPeakCapacityFactorMaximum, ...
    IsraelCombinedMonthlyPeakCapacityFactorAverage, IsraelCombinedMonthlyPeakCapacityFactorStd, ...
    resultsDirectory, nanmean(CFIsrael),IsraelWindDataYearVec);
