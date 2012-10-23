function optimizeMe = optimizeStations(stationWeightsConcatenated)
%optimizeStations optimizez best station choice, and station weights
%(amount of installed power in station A vs. station B etc.)
%   The optimization is done with fminsearch for a fixed number of
%   stations, for either 'all', 'winter', or 'summer' best fit to
%   electricity loads.
%   Initial guess is based on the rankStationCorrelation best guess, taking
%   a specific number of stations from there
% input vector is statioVec and weightVec concatenated

global input
global report
global optimizeFlag
global optimizeTheStations
% getting arrays out of struct 'input'
windEnergyNormalized = input.windEnergyNormalized;
tIsrael = input.tIsrael;
timeYearly = input.timeYearly;
Prated = input.Prated;
electricityNormalizedYearly = input.electricityNormalizedYearly;
minNan = input.minNan;
sortVector = input.sortVector;
stationVecOrig = input.stationVecOrig;

% choosing stations
N = length(stationWeightsConcatenated)/2;
if optimizeTheStations
    stationVec = round(stationWeightsConcatenated(1:N)); % fminsearch gives any float number
    stationVec(stationVec>31) = 31; % fminsearch has no limit to station number
else
    stationVec = stationVecOrig;
end
weights = stationWeightsConcatenated(N+1:end);

% normalizing weights to 1  
weights = weights / sum(weights);

% weighting them
for i=1:length(stationVec)
        weightedStations(i,:) = weights(i)*windEnergyNormalized(i,:);
end

% averaging them
[tVec,~,~,WindEnergyNormalizedAveraged, endVec, tfVec] = averageVectors(weightedStations,tIsrael,0);

% calculating metric - 
% fit to 'all', 'summer' or 'winter' electricity peak consumption
% TODO - ADD SOLAR!
[y, ~, ~, ~, ~, ~] = datevec(tVec);
IsraelWindDataYearVec = unique(y);
yearCounter = zeros(12,1);

% allocations
IsraelWindYearlyPeakCapacityFactorMinimum = NaN*ones(1,length(IsraelWindDataYearVec));
IsraelWindYearlyPeakCapacityFactorMaximum = NaN*zeros(1,length(IsraelWindDataYearVec));
IsraelWindYearlyPeakCapacityFactorAverage = zeros(1,length(IsraelWindDataYearVec));
IsraelWindYearlyPeakCapacityFactorStd = zeros(1,length(IsraelWindDataYearVec));
IsraelWindMonthlyPeakCapacityFactorMinimum = NaN*ones(1,12);
IsraelWindMonthlyPeakCapacityFactorMaximum = NaN*zeros(1,12);
IsraelWindMonthlyPeakCapacityFactorAverage = zeros(1,12);
IsraelWindMonthlyPeakCapacityFactorStd = zeros(1,12);

% yearly calculation
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
    
end

% debug
% plotting statistics of daily peak wind energy capacity
disp(['stations used: ' sprintf('%2.0d ',stationVec)])
disp(['weights used: ' sprintf('%2.2f ',weights)])
if report
    plotSimulatedPeakYear( maxDailyElectricity, IsraelWindDailyPeakCapacityFactor, CFIsrael(year-IsraelWindDataYearVec(1)+1), 'Israel combined', year, stationVec, weights)
end
% calculating and returning metric
switch optimizeFlag
    case 'all'
        optimizeMe = -mean(IsraelWindDailyPeakCapacityFactor);
    case 'winter'
        loc = [1:70,335:365];
        optimizeMe = -min(IsraelWindDailyPeakCapacityFactor(loc));
    case 'summer'
        loc = [150:340];
        optimizeMe = -min(IsraelWindDailyPeakCapacityFactor(loc));
end
disp(['daily ', optimizeFlag,' CC = ', num2str(-optimizeMe), '%'])
end

