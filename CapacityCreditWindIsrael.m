
function [stationVec] = CapacityCreditWindIsrael(stationVec, Ulimit, normPower, plotAll, totalStationMatFlag, optimize, optimizeFlag, optimizeTheStations,inputWeights)
% calculates the Capacity credit for wind, solar and combined solar+wind.
% this is done by the load method:
% 1. A average load year is calculated from the 20 year forecast data available from MEW through the IEC
% 2. The load data is ordered from high to low, and the sorting vector is saved.
% 3. The wind stations are chosen acording to combined minimum average wind speed and location based criterions
% 4. Each wind station is normalized to Katif data to give a "normPower" watt/m^2 energy density
% 5. Each station is multiplied by a representative power curve of a big wind turbine of 1 Mw rated power
% 6. For each station The capacity factor for each percent of loads is calculaed by P_avg/P_rated and plotted
% 7. The average capacity credit of all stations is summed to give an average station
% 8. The same is done for solar energy (only sde boker data)
% 9. The same is done for combined solar and wind (same rated power of both plants)

% TODOs:
% 1. add scenario titles to figures
% 2. add scenario numbers to fig and png file names (for creating movie?)
% 3. add selecting stations according to minimum years of data (or year
% span)
% 4. some kind of summer or winter index for instance - for high capacity
% credit at that time. Shavey Zion for instance rocks at the summer
%           Variable initialization
%-------------------------------------------------
global windDataDirectory 
global resultsDirectory 
optimizeFlagOld = optimizeFlag;
global optimizeFlag 
global optimizeTheStations
warning off
global debug; global report; global reportName; global joker
debug = 0; report = 0; joker = 0;
optimizeOrig = optimize;
reportName = ['results/report_' datestr(now) '.ps'];
if isunix
    windDataDirectory = '';
    resultsDirectory = [windDataDirectory 'results/'];
    electricityDataPath = '/home/hanan/Dropbox/MyPHDArticles/InProgress/IsraelDiurnalWind/YoelCohenData/';
end
if ispc
    windDataDirectory = 'D:\WindData\';
    resultsDirectory = [windDataDirectory 'results\'];
    electricityDataPath = 'C:\Users\hananel\Dropbox\MyPHDArticles\InProgress\IsraelDiurnalWind\YoelCohenData\';
end
if not(isdir(resultsDirectory))
    mkdir(resultsDirectory);
end

if nargin<9
    useBestCorrelatedNew=0;
else
    useBestCorrelatedNew=1;
end
months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
windDataYearVec = 2006:2011;
timeYearly = 1/6:1/6:365*24; % [hour] starts from 10 minutes - as the dates in the IMS data set so
% reading a typical big wind turbine power curve and normalizing the power.
% power curve chosen is Vestas 90 meter diameter 3 Mw wind turbine
[V,P,~] = VestasV903Mw102db();
P = P/max(P); Prated = 1;
% loading station meta-data
metaAll = loadMeta;
minNan = 6; %days with nans within a month that are tolerated for peak production statistics
%width = 0; % [hr] width of time around maximum consumption time to average the wind/solar power
percentVec = [0.001902587519026 0.1 0.5 1:100];
% 10 minute yearly peak is equivalent to
% 1/(24*365*6)*100=0.001902587519026% - kind of redundent. had to check :)
loc1 = find(percentVec==1,1);
loc30 = find(percentVec==30,1);
loc50 = find(percentVec==50,1);
loc100 = length(percentVec);
%           memory allocation
%-------------------------------------------------
EWYIsrael = double(zeros(length(windDataYearVec),length(timeYearly)));
weight = zeros(length(windDataYearVec),length(timeYearly));
CFper = NaN * ones(83,length(windDataYearVec),length(percentVec));
yearlyPeakCapacityFactorMinimum = zeros(83,6);
yearlyPeakCapacityFactorMaximum = zeros(83,6);
yearlyPeakCapacityFactorAverage = zeros(83,6);
yearlyPeakCapacityFactorStd = zeros(83,6);
dailyPeakCapacityFactor = zeros(365,1);
peakDailyTime = zeros(365,1);
stationWindEnergyNormalized = NaN * zeros(6*24*6*365,83);
stationWindEnergyNormalizedYearly = NaN * zeros(24*6*365,6);
monthlyPeakCapacityFactorMinimum = NaN*ones(83,12);
monthlyPeakCapacityFactorMaximum = NaN*zeros(83,12);
monthlyPeakCapacityFactorAverage = zeros(83,12);
monthlyPeakCapacityFactorStd = zeros(83,12);
IsraelWindMonthlyPeakCapacityFactorMinimum = NaN*ones(1,12);
IsraelWindMonthlyPeakCapacityFactorMaximum = NaN*zeros(1,12);
IsraelWindMonthlyPeakCapacityFactorAverage = zeros(1,12);
IsraelWindMonthlyPeakCapacityFactorStd = zeros(1,12);
IsraelSolarMonthlyPeakCapacityFactorMinimum = NaN*ones(1,12);
IsraelSolarMonthlyPeakCapacityFactorMaximum = NaN*zeros(1,12);
IsraelSolarMonthlyPeakCapacityFactorAverage = zeros(1,12);
IsraelSolarMonthlyPeakCapacityFactorStd = zeros(1,12);
IsraelCombinedMonthlyPeakCapacityFactorMinimum = NaN*ones(1,12);
IsraelCombinedMonthlyPeakCapacityFactorMaximum = NaN*zeros(1,12);
IsraelCombinedMonthlyPeakCapacityFactorAverage = zeros(1,12);
IsraelCombinedMonthlyPeakCapacityFactorStd = zeros(1,12);
CFtempMonthly = NaN*ones(83,12);
CF50Monthly = NaN*ones(83,12);
CF30Monthly = NaN*ones(83,12);
CF1Monthly = NaN*ones(83,12);
CFperMonthly = NaN*ones(83,6,length(percentVec));
Cor = zeros(83,6);
CorDaily = zeros(83,6,365);

%           matlab related commands
%-------------------------------------------------
more off;
% close all;
set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',12);
sp1 = 2; sp2 = 3; % subplot size for CC plots

%           0. if flagged for loading existing totalStationMat
%-------------------------------------------------
if totalStationMatFlag
    disp('Loading totalStationMat')
    load totalStationMat
    useBestCorrelated = useBestCorrelatedNew;
else
    
    %           1. calculate average load year
    %-------------------------------------------------
    [electricityNormalizedYearly,tYearly,electricityNormalizedDaily] = averageIECyear(0); %TODO - wierd number. 366 days?? ; to rerun - use averageIECyear(electricityDataPath)
    % interpolating to a 10 minute time vector, for 365 days.
    electricityNormalizedYearly = interp1(tYearly,electricityNormalizedYearly,timeYearly);
    
    %           2. ordering the load data from high to low
    %-------------------------------------------------
    [electricityNormalizedYearlySorted,sortVector]=sort(electricityNormalizedYearly,'descend');
    electricityNormalizedYearlySortedMonthly = zeros(12,31*24*6);
    sortVectorMonthly = zeros(12,31*24*6);
    % sorting per month
    for month=1:12
        dStart = datenum(0,month,0,0,0,0)*24*6+1;
        dEnd = min(datenum(0,month+1,0,0,0,0),365)*24*6;
        [electricityNormalizedYearlySortedMonthly(month,1:(dEnd-dStart+1)),sortVectorMonthly(month,1:(dEnd-dStart+1))] = ...
            sort(electricityNormalizedYearly(dStart:dEnd),'descend');
        sortVectorMonthly(month,1:(dEnd-dStart+1)) = sortVectorMonthly(month,1:(dEnd-dStart+1)) + dStart - 1;
        % plotting chronological method - month variation
        plotChronologicalMethod( timeYearly, electricityNormalizedYearly, sortVectorMonthly, ...
            sortVector, month, dStart, dEnd, 0.3, 0)
    end
    % plotting chronological method - year and day variation
    plotChronologicalMethod( timeYearly, electricityNormalizedYearly, sortVectorMonthly, ...
        sortVector, month, dStart, dEnd, 0.3, 1)
    if plotAll
        plotIEC(timeYearly,electricityNormalizedYearly,electricityNormalizedYearlySorted, sortVector, resultsDirectory);
    end
    
    %           3. choosing wind stations
    %-------------------------------------------------
    [stationVec, colOld] = filterStations(stationVec, Ulimit);
    stationNum1 = 0;
    stationVecOrig = stationVec;
    
    %              running over stations
    %-------------------------------------------------------------
    % additional allocations
    windEnergyNormalized = zeros(length(stationVecOrig),24*365*6*6);
    tIsrael = zeros(length(stationVecOrig),24*365*6*6);
    for Num=stationVecOrig
        
        % load data (windDataDirectory get's run over by load(matFile), also stationNum, who knows what else...)
        % windDataDirectory = '/home/hanan/Documents/measurements/';
        if isunix
            pathname = ['IMS-mat-files/',num2str(Num),'-',metaAll(Num).name,'/'];
        else
            pathname = ['IMS-mat-files\',num2str(Num),'-',metaAll(Num).name,'\'];
        end
        matFile = [pathname, 'Data_',num2str(Num),'.mat'];
        if exist(matFile)
            stationNum1 = stationNum1 + 1;
            disp(['loading ', num2str(stationNum1), ' ', metaAll(Num).name , ' matFile'])
            load(matFile);
            if isunix
                resultsDirectory = [windDataDirectory 'results/'];
            else
                resultsDirectory = [windDataDirectory 'results\'];
            end
            [y, m, d, h, mi] = datevec(t);
            stationNum = stationNum1;
            
            % 4. Normalizing power of stations to 'normPower' watt/m^2 (Katif=390 W/m^2 @ 60 m)
            normFactor = 0.5*1.2*nansum(U.^3)/length(U); % W/m^2
            Unormalized(1:length(U),Num) = U * (normPower/normFactor)^(1/3);
            
            % 4* - summing daily data for visualization purposes
            UDailyNormalized(:,:,stationNum1) = UDaily * (normPower/normFactor)^(1/3);
            
            % 5. multiplying by Prated Mw power curve
            stationWindEnergyNormalized(1:length(U),Num) = interp1(V,P,Unormalized(1:length(U),Num));
            
            % 5* - multiplying the daily data
            stationWindEnergyNormalizedDaily(:,:,stationNum1) = interp1(V,P,UDailyNormalized(:,:,stationNum1));
            
            % 6. the Capacity credit is calculated, for each year of data
            yearVec=unique(y);
            if length(yearVec)>7
                [y, m, d, h, mi, s] = datevec(t);
                yearVec=[unique(y)]';
            end
            yearVec(yearVec==2012) = [];
            
            %clear CFtemp CF50 CF30 Cor
            CFtemp = NaN*ones(1,length(windDataYearVec));
            CF50 = NaN*ones(1,length(windDataYearVec));
            CF30 = NaN*ones(1,length(windDataYearVec));
            CF1 = NaN*ones(1,length(windDataYearVec));
            Cor = NaN*ones(1,length(windDataYearVec));
            
            % point 7. summing all stations - for representing all Israel. later on looking for the overlapping years
            windEnergyNormalized(stationNum,1:length(stationWindEnergyNormalized)) = stationWindEnergyNormalized(:,Num);
            tIsrael(stationNum,1:length(t)) = t;
            
            % running over years of data for each station
            % ----------------------------------------------------------------
            yearCounter = zeros(12,1); yearCounter1 = zeros(12,1); yearCounter30 = zeros(12,1);
            yearCounter50 = zeros(12,1); yearCounter100 = zeros(12,1);
            % need "ismatlab" switch - for transposing yearVec for octave
            for year=yearVec
                % taking yearly data
                temp = stationWindEnergyNormalized(y==year,Num);
                stationWindEnergyNormalizedYearly(1:length(temp),year-windDataYearVec(1)+1) = temp;
                UYearly = U(y==year);
                % interpolating to timeYearly
                tOneYearWind = t(find(y==year));
                tOneYearWind = tOneYearWind-tOneYearWind(1); %[day]
                % elimiating repeated t
                [b,bi,bj] = unique(tOneYearWind);
                stationWindEnergyNormalizedYearlyInterp = interp1(b', stationWindEnergyNormalizedYearly(bi, year - windDataYearVec(1) + 1), timeYearly / 24);
                UYearlyInterp = interp1(b',UYearly(bi),timeYearly/24);
                CF(Num,year-windDataYearVec(1)+1) = nansum(stationWindEnergyNormalizedYearlyInterp)/length(stationWindEnergyNormalizedYearlyInterp)/Prated; % missing multiplication by deltaT. works because all is normalized
                
                % arranging the wind data in the sorted order of the load
                stationWindEnergyNormalizedYearlyInterpSort = stationWindEnergyNormalizedYearlyInterp(sortVector);
                % comparing for 1%-100% top loads
                percentCounter = 0;
                for percent=percentVec
                    percentCounter = percentCounter+1;
                    loc = round(24*365*percent/100*6); % the 6 factor is for the 1 hour resolution for the load data vs. the 10 minute wind data
                    CFper(Num,year-windDataYearVec(1)+1,percentCounter) = nansum(stationWindEnergyNormalizedYearlyInterpSort(1:loc))/loc/Prated;
                end
                CFtemp(year-windDataYearVec(1)+1) = CF(Num,year-windDataYearVec(1)+1);
                CF50(year-windDataYearVec(1)+1) = CFper(Num,year-windDataYearVec(1)+1,loc50);
                CF30(year-windDataYearVec(1)+1) = CFper(Num,year-windDataYearVec(1)+1,loc30);
                CF1(year-windDataYearVec(1)+1) = CFper(Num,year-windDataYearVec(1)+1,loc1);
                
                % monthly Capacity Credit
                for month=1:12
                    % arranging the wind data in the sorted order of the load
                    dStart = datenum(0,month,0,0,0,0)+1;
                    dEnd = min(datenum(0,month+1,0,0,0,0),365);
                    lMonth = (dEnd-dStart)*24*6; % verify
                    stationWindEnergyNormalizedYearlyInterpSortMonthly = stationWindEnergyNormalizedYearlyInterp(sortVectorMonthly(month,1:lMonth));
                    % comparing for 1%-100% top loads
                    nan1=0; nan30=0; nan50=0; nan100=0;
                    percentCounter = 0;
                    for percent=percentVec
                        percentCounter = percentCounter+1;
                        loc = round(lMonth*percent/100); % the 6 factor is for the 1 hour resolution for the load data vs. the 10 minute wind data
                        % checking if there are too many nans in the monthly
                        % data
                        nanNum = sum(isnan(stationWindEnergyNormalizedYearlyInterpSortMonthly(1:loc))); % not used at the moment
                        if nanNum/loc<=minNan/30
                            CFperMonthly(Num,year-windDataYearVec(1)+1,percentCounter) = nansum(stationWindEnergyNormalizedYearlyInterpSortMonthly(1:loc))/loc/Prated;
                            if percent==1
                                nan1 = 1;
                                yearCounter1(month) = yearCounter1(month)+1;
                            end
                            if percent==30
                                nan30 = 1;
                                yearCounter30(month) = yearCounter30(month)+1;
                            end
                            if percent==50
                                nan50 = 1;
                                yearCounter50(month) = yearCounter50(month)+1;
                            end
                            if percent==100
                                nan100 = 1;
                                yearCounter100(month) = yearCounter100(month)+1;
                            end
                        end
                        pp(percentCounter) = CFperMonthly(Num,year-windDataYearVec(1)+1,percentCounter);
                    end
                    if nan100
                        CFtempMonthly(stationNum,month) = nansum([CFtempMonthly(stationNum,month)*(yearCounter100(month)-1), ...
                            CFperMonthly(Num,year - windDataYearVec(1) + 1,loc100)])/yearCounter100(month); end
                    if nan50
                        CF50Monthly(stationNum,month)   = nansum([CF50Monthly(stationNum,month)*(yearCounter50(month)-1), ...
                            CFperMonthly(Num,year - windDataYearVec(1) + 1,loc50)])/yearCounter50(month); end
                    if nan30
                        CF30Monthly(stationNum,month)   = nansum([CF30Monthly(stationNum,month)*(yearCounter30(month)-1), ...
                            CFperMonthly(Num,year - windDataYearVec(1) + 1,loc30)])/yearCounter30(month); end
                    if nan1
                        CF1Monthly(stationNum,month)    = nansum([CF1Monthly(stationNum,month)*(yearCounter1(month)-1), ...
                            CFperMonthly(Num,year - windDataYearVec(1) + 1,loc1)])/yearCounter1(month); end
                    if debug
                        figure(131); subplot(212); plot(percentVec, pp); title([months(month), ' month count ', num2str(yearCounter50(month))])
                        subplot(211); plot(1:12,CF50Monthly(stationNum,:),'r',1:12,CF30Monthly(stationNum,:),'g');
                        title(metaAll(Num).name);
                    end
                end
                % calculating correlation
                % YOEL - we may want to test correlation per month - I think
                % going for the summer and transition months makes the most
                % sense. showing the case there
                % a. the minimum capcaity credit of all wind stations is higher
                % then zero
                % b. by designating Prated_per_station we can get combinations
                % that are even better
                % c. combination with solar increaces capacity credit at peak
                % times
                [Cor, CorDaily] = recordStationCorrelation( ...
                    timeYearly, stationWindEnergyNormalizedYearlyInterp, ...
                    electricityNormalizedYearly, year, windDataYearVec, Num, ...
                    CorDaily, Cor);
                
                % recording maximum,minimum,average and std of capacity factor
                % at peak consumption for each station
                % TODO - why dont i return and update yearCounter??
                [monthlyPeakCapacityFactorMinimum, monthlyPeakCapacityFactorMaximum,monthlyPeakCapacityFactorAverage, monthlyPeakCapacityFactorStd, ...
                    yearlyPeakCapacityFactorMinimum, yearlyPeakCapacityFactorMaximum, yearlyPeakCapacityFactorAverage, yearlyPeakCapacityFactorStd, ...
                    peakDailyTime, dailyPeakCapacityFactor, yearCounter] = ...
                    recordPeakProductionStatistics(timeYearly, electricityNormalizedYearly, stationWindEnergyNormalizedYearlyInterp, ...
                    Prated, yearCounter, windDataYearVec, Num, year, minNan, ...
                    monthlyPeakCapacityFactorMinimum, monthlyPeakCapacityFactorMaximum, ...
                    monthlyPeakCapacityFactorAverage, monthlyPeakCapacityFactorStd, ...
                    yearlyPeakCapacityFactorMinimum, yearlyPeakCapacityFactorMaximum, ...
                    yearlyPeakCapacityFactorAverage, yearlyPeakCapacityFactorStd, metaAll(Num).name);
                
                if plotAll % plotting station data, inside year loop
                    % time series of each year
                    tic
                    disp([' plotting timeline for ' num2str(year)])
                    plotTimeLine(metaAll(Num).name,Num,year-windDataYearVec(1)+1,year,CFtemp(year-windDataYearVec(1)+1),timeYearly,electricityNormalizedYearly, ...
                        stationWindEnergyNormalizedYearlyInterp,dailyPeakCapacityFactor,peakDailyTime,Prated,resultsDirectory);
                    toc
                    
                    % Capacity credit of each year vs. precent of maximum
                    % electricity usage
                    disp([' plotting capacity credit for ' num2str(year)])
                    plotCC(percentVec, sp1,sp2,Num,year,year-windDataYearVec(1)+1,CFper(Num,year-windDataYearVec(1)+1,:),1)
                end
                
                if debug
                    figure(8888+Num);
                    title(metaAll(Num).name);
                    subplot(211) ; hold on;
                    plot(timeYearly/24, moving(stationWindEnergyNormalizedYearlyInterp,30*6*24,'nanmean'))
                    subplot(212) ; hold on;
                    plot(timeYearly/24, moving(UYearlyInterp,30*6*24,'nanmean'))
                end
                
                if report
                    % making a representative 1 week plot for summer (and
                    % winter?)
                    % using station 8 - shavei zion, which shows a very nice
                    % distinct summer pattern (alternatives - 34 Gilboa, 13
                    % Eshar)
                    if and(Num==13,year==2011)
                        plotWeek( timeYearly, stationWindEnergyNormalizedYearlyInterp, electricityNormalizedYearly, ...
                            [num2str(Num),'-',metaAll(Num).name], year, 10001, 10002)
                    end
                end
            end
            
            if plotAll % plotting station data outside of year loop
                for year=windDataYearVec
                    if ~(yearVec==year)
                        subplot(sp1,sp2,year-windDataYearVec(1)+1);
                        axis off; box off;
                    end
                end
                suptitle([meta.name ' , capacity Credit for precentage of peak demand'])
                set(gcf,'Units','normalized')
                set(gcf,'Position',[0.126647 -0.294271 1.13031 1.10417])
                print([resultsDirectory 'CCfull_', strrep(meta.name,' ',''),'.png'],'-dpng');
                saveas(gcf, [resultsDirectory 'CCfull_', strrep(meta.name,' ','')], 'fig')
                
                % plotting statistics of daily peak wind energy capacity
                disp([' plotting peak statistics for ' num2str(year)])
                plotPeak(metaAll(Num).name, Num, yearlyPeakCapacityFactorMinimum, yearlyPeakCapacityFactorMaximum, ...
                    yearlyPeakCapacityFactorAverage, yearlyPeakCapacityFactorStd, ...
                    monthlyPeakCapacityFactorMinimum, monthlyPeakCapacityFactorMaximum, ...
                    monthlyPeakCapacityFactorAverage, monthlyPeakCapacityFactorStd, ...
                    resultsDirectory, nanmean(CFtemp),yearVec);
                
                % plotting yearly 30 and 50 markers, plus correlation
                plot3050CC(Num, meta.name ,resultsDirectory, windDataYearVec,CFtemp,CF50,CF30,CF1,Cor,yearVec)
                
                % plotting each months's daily profile (avg. over al years of data) compared to load
                plotMonthlyAverageNormalizedWindEnergy( Num, meta.name, electricityNormalizedDaily, UDaily, months, resultsDirectory )
                
            end
        end
    end
    
    % save debug mat
    save totalStationMat
end
optimizeFlag = optimizeFlagOld;
% 6.5 - calculate station ranking based on correlation so that I can try to
% cross correlt the stations themselves!
% basis - look for stations with most days that are non correlated (look at
% correlation of the electricity-correlation matrix
if optimizeFlag==0
    optimizeFlag = 'summer';
end
[bestStationSetNum, bestStationSetVec, sortedStationVec, bestStationSetWeight, finalRank, a, b] = rankStationCorrelation(CorDaily, Cor, stationVecOrig,optimizeFlag, metaAll);

% 7. average capacity credit for each year for all of Israel, including solar

% TODO -
% looking for overlapping years and averaging wind data and CC
% at the moment - taking best correlation can offer!
% for year = 2010, yearNum=5
yearNum = 5;
stations = sortedStationVec{yearNum}(bestStationSetVec{yearNum});

optimize = optimizeOrig;
if optimize
    % using fminsearch to get best optimization, with the narrowing of stations
    % based on the ranking system, and optimization of all, winter or summer
    % months
    % optimizationResults return stationVec and weightVec concatenated
    % starting with best 2 stations
    for stationNumOptimize = 2:10
        % generating input and initial guess
        % Iinitial guess
        weights{stationNumOptimize} = ones(1,stationNumOptimize)/stationNumOptimize;
        %weights{stationNumOptimize}(stationNumOptimize) = 1;
        % input
        x = [stations(1:stationNumOptimize) weights(1:stationNumOptimize)];
        global input
        input.windEnergyNormalized = windEnergyNormalized(stations(1:stationNumOptimize),:);
        input.tIsrael = tIsrael(stations(1:stationNumOptimize),:);
        input.timeYearly = timeYearly;
        input.Prated = Prated;
        input.electricityNormalizedYearly = electricityNormalizedYearly;
        input.minNan = minNan;
        input.sortVector = sortVector;
        input.stationVecOrig = stations(1:stationNumOptimize);
        Options = optimset('fmincon');
        Options.TolX = 1e-12;
        Options.TolFun = 1e-12;
        Options.UseParallel = 'always';
        Options.DiffMinChange = 1.0000e-02;
        joker = stationNumOptimize;
        optimizeWeights(weights{stationNumOptimize});
        report = 0;
        if optimizeTheStations
            [optimizationResults optimizeMe(stationNumOptimize)] = fminsearch(@(x) optimizeStations(x), x, Options);
        else
            % TODO - add the constraint that weights*ones(length(weights),1)=1 = 
            [X, FVAL] = fmincon(@(x) optimizeWeights(x),weights{stationNumOptimize}',[],[],ones(1,length(weights)),1,zeros(1,length(weights)),ones(1,length(weights)),[],Options);
            weights{stationNumOptimize} = X; 
            optimizeMe(stationNumOptimize) = FVAL;%  fminsearch(@(x) optimizeWeights(x), weights, Options);
        end
        if optimizeTheStations
            N = length(optimizationResults)/2;
            stationVec = optimizationResults(1:N);
            weights = optimizationResults(N+1:end);
        end
        for i=1:stationNumOptimize
            weightedStations(i,:) = weights{stationNumOptimize}(i)*windEnergyNormalized(stations(i),:);
        end
        report = 1; joker = 1+stationNumOptimize;
        optimizeWeights(weights{stationNumOptimize});
        report = 0; joker = 0;
        disp('saving')
        disp(weights)
        disp(stationVecOrig(stations(1:stationNumOptimize)))
        save lastStand optimizeMe weights stations stationNumOptimize stationVecOrig metaAll
    end
   
    [tVec,~,~,WindEnergyNormalizedAveraged, endVec, tfVec] = averageVectors(weightedStations,tIsrael(stations(1:stationNumOptimize),:),0);
    [tVecM,~,~,CFMonthlyIsrael, endVecM, monthsVec] = averageVectors(CFtempMonthly(stations(1:stationNumOptimize),:),1:12,0);
    [tVecM,~,~,CF50MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF50Monthly(stations(1:stationNumOptimize),:),1:12,0);
    [tVecM,~,~,CF30MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF30Monthly(stations(1:stationNumOptimize),:),1:12,0);
    [tVecM,~,~,CF1MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF1Monthly(stations(1:stationNumOptimize),:),1:12,0);
else
    if useBestCorrelated
         % no optimization at the moment - just same weight for each station
        load bestSummer.mat
        stationNum = length(weights);
        stationNumOptimize = stationNum;
        counter = 0;
        for i=stations
            counter = counter+1;
            weightedStations(counter,:) = weights(counter)*windEnergyNormalized((i),:);
        end
        [tVec,~,~,WindEnergyNormalizedAveraged, endVec, tfVec] = averageVectors(weightedStations,tIsrael(stations,:),0);
        % old version:
        % [tVec,~,~,WindEnergyNormalizedAveraged, endVec, tfVec] = averageVectors(windEnergyNormalized,tIsrael,0);
        [tVecM,~,~,CFMonthlyIsrael, endVecM, monthsVec] = averageVectors(CFtempMonthly(1:stationNum,:),1:12,0);
        [tVecM,~,~,CF50MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF50Monthly(1:stationNum,:),1:12,0);
        [tVecM,~,~,CF30MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF30Monthly(1:stationNum,:),1:12,0);
        [tVecM,~,~,CF1MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF1Monthly(1:stationNum,:),1:12,0);
        disp('took best optimization for summer - highest correlated 10 stations!')
        disp('best stations and weights --->')
        disp([stationVecOrig(stations)' , weights])
        plotStationLocations(stationVecOrig(stations(1:stationNum)),jet(83),weights)

    else
        % no optimization at the moment - just same weight for each station
        weights = ones(length(stationVecOrig),1)/length(stationVecOrig); %bestStationSetWeight{yearNum}';
        for i=1:length(stationVecOrig)
            weightedStations(i,:) = weights(i)*windEnergyNormalized(i,:);
        end
        [tVec,~,~,WindEnergyNormalizedAveraged, endVec, tfVec] = averageVectors(weightedStations,tIsrael,0);
        % old version:
        % [tVec,~,~,WindEnergyNormalizedAveraged, endVec, tfVec] = averageVectors(windEnergyNormalized,tIsrael,0);
        [tVecM,~,~,CFMonthlyIsrael, endVecM, monthsVec] = averageVectors(CFtempMonthly(1:stationNum,:),1:12,0);
        [tVecM,~,~,CF50MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF50Monthly(1:stationNum,:),1:12,0);
        [tVecM,~,~,CF30MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF30Monthly(1:stationNum,:),1:12,0);
        [tVecM,~,~,CF1MonthlyIsrael, endVecM, monthsVec] = averageVectors(CF1Monthly(1:stationNum,:),1:12,0);
    end
end  

if report
    plotTimeLines(windEnergyNormalized,tIsrael,endVec, tfVec,colOld, WindEnergyNormalizedAveraged,tVec, stationVecOrig);
end

% loading solar data
[NIPDailyNorm_summer,NIPDailyNorm_winter,tDailySolar,axSolar,NIP,NIPDaily,dSolar,mSolar,ySolar,hSolar,miSolar] = loadIMSradiation(75,5003,windDataDirectory,0);
% add solar plot to week typical plot
%TODO plotWeek

% some more allocations
IsraelSolarDataYearVec = unique(ySolar);
[y, m, d, h, minute, sec] = datevec(tVec);
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
yearCounterSolar = zeros(12,1);
yearCounterCombined = zeros(12,1);

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
        peakDailyTime, IsraelWindDailyPeakCapacityFactor, maxDailyElectricity, yearCounter] = ...
        recordPeakProductionStatistics(timeYearly, electricityNormalizedYearly, WindEnergyNormalizedAveragedYearlyInterp, ...
        Prated, yearCounter, IsraelWindDataYearVec, 1, year, minNan, ...
        IsraelWindMonthlyPeakCapacityFactorMinimum, IsraelWindMonthlyPeakCapacityFactorMaximum, ...
        IsraelWindMonthlyPeakCapacityFactorAverage, IsraelWindMonthlyPeakCapacityFactorStd, ...
        IsraelWindYearlyPeakCapacityFactorMinimum, IsraelWindYearlyPeakCapacityFactorMaximum, ...
        IsraelWindYearlyPeakCapacityFactorAverage, IsraelWindYearlyPeakCapacityFactorStd, 'Israel wind');
    % simulation of 15000Mw/1500Mw scenario
    joker = 5;
    plotSimulatedPeakYear( maxDailyElectricity, IsraelWindDailyPeakCapacityFactor, CFIsrael(year-IsraelWindDataYearVec(1)+1), 'Israel combined', year)
    if year==2011
        plotWeek( timeYearly, WindEnergyNormalizedAveragedYearlyInterp, electricityNormalizedYearly...
            , 'Israel combined', year, 10003, 10004)
    end
    % Solar
    [IsraelSolarMonthlyPeakCapacityFactorMinimum, IsraelSolarMonthlyPeakCapacityFactorMaximum,IsraelSolarMonthlyPeakCapacityFactorAverage, IsraelSolarMonthlyPeakCapacityFactorStd, ...
        IsraelSolarYearlyPeakCapacityFactorMinimum,  IsraelSolarYearlyPeakCapacityFactorMaximum, IsraelSolarYearlyPeakCapacityFactorAverage,  IsraelSolarYearlyPeakCapacityFactorStd, ...
        peakDailyTime, IsraelSolarDailyPeakCapacityFactor, yearCounterSolar] = ...
        recordPeakProductionStatistics(timeYearly, electricityNormalizedYearly, NIPnormalizedYearlyInterp, ...
        Prated, yearCounterSolar, IsraelSolarDataYearVec, 1, year, minNan, ...
        IsraelSolarMonthlyPeakCapacityFactorMinimum, IsraelSolarMonthlyPeakCapacityFactorMaximum, ...
        IsraelSolarMonthlyPeakCapacityFactorAverage, IsraelSolarMonthlyPeakCapacityFactorStd, ...
        IsraelSolarYearlyPeakCapacityFactorMinimum, IsraelSolarYearlyPeakCapacityFactorMaximum, ...
        IsraelSolarYearlyPeakCapacityFactorAverage, IsraelSolarYearlyPeakCapacityFactorStd, 'Israel solar');
    
    % Combined
    [IsraelCombinedMonthlyPeakCapacityFactorMinimum, IsraelCombinedMonthlyPeakCapacityFactorMaximum,IsraelCombinedMonthlyPeakCapacityFactorAverage, IsraelCombinedMonthlyPeakCapacityFactorStd, ...
        IsraelCombinedYearlyPeakCapacityFactorMinimum,  IsraelCombinedYearlyPeakCapacityFactorMaximum, IsraelCombinedYearlyPeakCapacityFactorAverage,  IsraelCombinedYearlyPeakCapacityFactorStd, ...
        peakDailyTime, IsraelCombinedDailyPeakCapacityFactor, yearCounterCombined] = ...
        recordPeakProductionStatistics(timeYearly, electricityNormalizedYearly, WindEnergyNIPnormalizedYearlyInterp, ...
        Prated, yearCounterCombined, IsraelWindDataYearVec, 1, year, minNan, ...
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
    plotCC(percentVec, sp1,sp2,Num,year,year-IsraelWindDataYearVec(1)+1, ...
        CFperIsrael(year-IsraelWindDataYearVec(1)+1,:),0, CFperIsraelS(year-IsraelWindDataYearVec(1)+1,:), CFperIsraelSW(year-IsraelWindDataYearVec(1)+1,:));
end
suptitle('Israel scenario - capacity Credit for precentage of peak demand')
set(gcf,'Units','normalized')
set(gcf,'Position',[0.126647 -0.294271 1.13031 1.10417])
legend('wind','solar','wind+solar','Location','Best')
saveas(gcf, [resultsDirectory 'IsraelCC_'], 'fig')

% plot spectacular presentation front image
plotDailyWindSolarElectricityVisualization(electricityNormalizedDaily,stationWindEnergyNormalizedDaily,NIPDaily)

% plotting combined stations

figure(400)
bar((ones(3,1)*IsraelWindDataYearVec)',[CFIsrael50;CFIsraelS50;CFIsraelSW50]'); hold on;
legend('Wind','Solar','Wind+Solar','Location','North');
minCF50 = floor(min(min([CFIsrael50;CFIsraelS50;CFIsraelSW50]))*10)/10-0.03; maxCF50 = ceil(max(max([CFIsrael50;CFIsraelS50;CFIsraelSW50]))*100)/100+0.01;
[h_ax,h_l] = plot2ylinscales([min(IsraelWindDataYearVec),max(IsraelWindDataYearVec)],[minCF50 maxCF50],1/nanmean(CFIsrael));
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
bar((ones(3,1)*IsraelWindDataYearVec)',[CFIsrael30;CFIsraelS30;CFIsraelSW30]'); hold on;
legend('Wind','Solar','Wind+Solar','Location','North');
minCF30 = floor(min(min([CFIsrael30;CFIsraelS30;CFIsraelSW30]))*10)/10-0.03; maxCF30 = ceil(max(max([CFIsrael30;CFIsraelS30;CFIsraelSW30]))*100)/100+0.01;
[h_ax,h_l] = plot2ylinscales([min(IsraelWindDataYearVec),max(IsraelWindDataYearVec)],[minCF30 maxCF30],1/nanmean(CFIsrael));
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
bar((ones(3,1)*IsraelWindDataYearVec)',[CFIsrael1;CFIsraelS1;CFIsraelSW1]'); hold on;
legend('Wind','Solar','Wind+Solar','Location','North');
minCF1 = floor(min(min([CFIsrael1;CFIsraelS1;CFIsraelSW1]))*10)/10-0.03; maxCF1 = ceil(max(max([CFIsrael1;CFIsraelS1;CFIsraelSW1]))*100)/100+0.01;
[h_ax,h_l] = plot2ylinscales([min(IsraelWindDataYearVec),max(IsraelWindDataYearVec)],[minCF1 maxCF1],1/nanmean(CFIsrael));
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
plotMonthlyCC('Israel wind', CFMonthlyIsrael, CF50MonthlyIsrael, CF30MonthlyIsrael, CF1MonthlyIsrael, nanmean(CFIsrael))

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

if report % generate .ps report file
    figure;
    text(0,0.5,{['Report generated on the ' datestr(now)],...
        '-----------------------------------------', ...
        ['Minimum interanuual measured wind speed = ' num2str(Ulimit) ' m/s'], ...
        ['Power density = ' num2str(normPower) ' W/m^2']},'fontsize',18);
    axis([-0.1 0.7 0.4 0.6])
    axis off
    print('-dpsc','-append',reportName)
    figure(1); % stations locations
    print('-dpsc','-append',reportName)
    figure(3000);
    print('-dpsc','-append',reportName)
    
    figure(997); % chronological method - yearly
    print('-dpsc','-append',reportName)
    figure(998); % chronological method - monthly
    print('-dpsc','-append',reportName)
    figure(999); % "chronological method" - daily peak
    print('-dpsc','-append',reportName)
    
    figure;
    text(0,0.5,{'Typical temporal trend lines',...
        '-----------------------------------------'},'fontsize',18);
    axis([-0.1 0.7 0.4 0.6])
    axis off
    print('-dpsc','-append',reportName)
    figure(10001); % typical summer week
    print('-dpsc','-append',reportName)
    figure(10002); % typical winter week
    print('-dpsc','-append',reportName)
    figure(10003); % typical summer week
    print('-dpsc','-append',reportName)
    figure(10004); % typical winter week
    print('-dpsc','-append',reportName)
    
    figure;
    text(0,0.5,{'Simple simulation analysis results',...
        '-----------------------------------------'},'fontsize',18);
    axis([-0.1 0.7 0.4 0.6])
    axis off
    print('-dpsc','-append',reportName)
    figure(10015); % yearly total Israel chronological method - wind, solar, combined
    print('-dpsc','-append',reportName)
    figure(400); % bar plot - year vs. 30%CC - wind, solar, combined
    print('-dpsc','-append',reportName)
    figure(401); % bar plot - year vs. 50%CC - wind, solar, combined
    print('-dpsc','-append',reportName)
    figure(2490); % Israel - inter-annual monthly averaged daily peak of capacity credit
    print('-dpsc','-append',reportName)
    figure(10008); % Israel - inter-annual monthly and yearly capacity credit
    print('-dpsc','-append',reportName)
    for year=IsraelWindDataYearVec
        figure(10010+year+joker); %simulation of yearly scenario
        print('-dpsc','-append',reportName)
    end
    
    figure;
    text(0,0.5,{'Optimization analysis results',...
        '-----------------------------------------'},'fontsize',18);
    axis([-0.1 0.7 0.4 0.6])
    axis off
end

% best summer optimization - with old IEC average years (averag over all
% 2011-2030 forecast years)
% station num | station-num-from 1:83 list | weights
%   18.0000          7.0000            0.1229
%    13.0000          5.0000            0.2553
%    83.0000         26.0000            0.2634
%    16.0000          6.0000            0.0664
%    53.0000         16.0000            0.0729
%    72.0000          22.0000            0.1233
%     4.0000           2.0000            0.0158
%    11.0000           3.0000            0.0797
%    12.0000           4.0000            0.0004
%    71.0000          21.0000            0.0000
