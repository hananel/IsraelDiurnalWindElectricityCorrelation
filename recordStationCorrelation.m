function [Cor, CorDaily] = recordStationCorrelation( ...
                           timeYearly, stationWindEnergyNormalizedYearlyInterp, ...
                           electricityNormalizedYearly, year, windDataYearVec, Num, ...
                           CorDaily, Cor)
% rankStationCorrelation will rank each stations correlations to electricity
% production. Throughout the year, daily correlations, and calculate a "rank" 
% for the summer and winter

% yearly correlation
stationWindEnergyNormalizedYearlyInterp_cor = stationWindEnergyNormalizedYearlyInterp;
stationWindEnergyNormalizedYearlyInterp_cor(isnan(stationWindEnergyNormalizedYearlyInterp_cor)) = 0; % corr doesn't work with NaNs
Cor(Num, year-windDataYearVec(1)+1) = corr(stationWindEnergyNormalizedYearlyInterp_cor',electricityNormalizedYearly');

% daily correlation
for i=0:364
    loc = and(timeYearly/24>i,timeYearly/24<i+1);
    stationWindEnergyNormalizedDailyInterp = stationWindEnergyNormalizedYearlyInterp(loc);
    electricityNormalizedDaily = electricityNormalizedYearly(loc);
    stationWindEnergyNormalizedDailyInterp(isnan(stationWindEnergyNormalizedDailyInterp)) = 0; % corr doesn't work with NaNs
    CorDaily(Num, year-windDataYearVec(1)+1, i+1) = corr(stationWindEnergyNormalizedDailyInterp',electricityNormalizedDaily');
end
end

