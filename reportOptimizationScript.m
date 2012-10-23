% script for analyzing optimization results, and producing a report
warning off
close all

% setting
doWinter = 0;
doSummer = 1;
useOnlyLastResult = 1;
timeYearly = 1/6:1/6:365*24; % [hour] starts from 10 minutes - as the dates in the IMS data set so

%the easiest thing would have been to call CapacityCredit...m with the
%stations and weights for all optimization solutions ...
    %           1. calculate average load year
    %-------------------------------------------------
    [electricityNormalizedYearly,tYearly,electricityNormalizedDaily] = averageIECyear(0); %TODO - wierd number. 366 days?? ; to rerun - use averageIECyear(electricityDataPath)
    % interpolating to a 10 minute time vector, for 365 days.
    electricityNormalizedYearly = interp1(tYearly,electricityNormalizedYearly,timeYearly);
    
    
if doWinter
    % winter analysis
    disp('loading  winter results')
    load totalStationMat_optimization_winter_2to16.mat
    load optimizationWeight_winter_2to16
    stationNumOptimizeVec = 1:16;
    
    if useOnlyLastResult
        stationNumOptimizeVec = stationNumOptimizeVec(end);
    end
    
    % running through optimization cases
    for stationNumOptimize=stationNumOptimizeVec
        simulationName = ['Israel optimized winter #',num2str(stationNumOptimize)];
        colors = jet(stationNumOptimize);
        
        % use optimization based stations
        clear weightedStations
        fig1 = figure;
        fig2 = figure;
        for i=1:stationNumOptimize
            if stationNumOptimize==1
                weightedStations = windEnergyNormalized(stations(1),:);
                weights{1} = 1;
            else
                weightedStations(i,:) = weights{stationNumOptimize}(i)*windEnergyNormalized(stations(i),:);
            end
            if find(tIsrael(stations(i),:)==0,1)
                loc = find(tIsrael(stations(i),:)==0,1);
                figure(fig1); hold on; plot(tIsrael(stations(i),1:loc-1),windEnergyNormalized(stations(i),1:loc-1),'.','color',colors(i,:))
                figure(fig2); hold on; plot(tIsrael(stations(i),1:loc-1),weightedStations(i,1:loc-1),'.','color',colors(i,:))
            else
                figure(fig1); hold on; plot(tIsrael(stations(i),:),windEnergyNormalized(stations(i),:),'.','color',colors(i,:))
                figure(fig2); hold on; plot(tIsrael(stations(i),:),weightedStations(i,:),'.','color',colors(i,:))
            end
        end
        [tVec,~,~,WindEnergyNormalizedAveraged, endVec, tfVec] = averageVectors(weightedStations,tIsrael(stations(1:stationNumOptimize),:),0);
        figure(fig1); hold on; plot(tVec,WindEnergyNormalizedAveraged,'k');
        datetick
        figure(fig1); legend({metaAll(stationVecOrig(stations(1:stationNumOptimize))).name,'averaged'})
        figure(fig2); hold on; plot(tVec,WindEnergyNormalizedAveraged,'k');
        datetick
        figure(fig2); legend({metaAll(stationVecOrig(stations(1:stationNumOptimize))).name,'averaged'})
        
        % re-run analysis for the weights{i} and stationVecOrig(stations(1:i))
        optimizationReRunAnalysisScript
        
        % plot map of stations and weights - as a circle of proportional size
        % to weights{i}(stationNum)
        figure(9);
        plotStationLocations(stationVecOrig(stations(1:stationNumOptimize)),jet(83),weights{stationNumOptimize})
        print('-dpsc','-append',reportName)
        
        % calculate optimization parameter again
        if stationNumOptimize>1
            global input
            input.windEnergyNormalized = windEnergyNormalized(stations(1:stationNumOptimize),:);
            input.tIsrael = tIsrael(stations(1:stationNumOptimize),:);
            input.timeYearly = timeYearly;
            input.Prated = Prated;
            input.electricityNormalizedYearly = electricityNormalizedYearly;
            input.minNan = minNan;
            input.sortVector = sortVector;
            input.stationVecOrig = stations(1:stationNumOptimize);
            optimizationNormalized(stationNumOptimize) = -100*optimizeWeights(weights{stationNumOptimize})/CFIsrael(end-1); % [%] - reference year is 2011
        end
        
        % script for saving plots to report as .ps file
        reportDir = ['~/Documents/diurnalIsrael/results/reportWinter',num2str(stationNumOptimize),'/'];
        reportName = ['optimizationReport_winter_2to16_' num2str(stationNumOptimize) '.ps'];
        mkdir(reportDir)
        printOptimizationReport
        close all
        
    end
    
    % optimization convergence
    figure(1); hold on;
    plot(2:16,optimizationNormalized(2:16),'*-')
    title({'Minimum winter daily-peak capacity credit','as a function of stations used in optimization'})
    xlabel('Number of top ranking stations used')
    ylabel('Minimum winter daily-peak CC/CF [%]')
    axis tight
    ticks = get(gca,'YTick');
    for i=1:length(ticks)
        leg{i} = [num2str(ticks(i)),'%'];
    end
    set(gca,'YTickLabel',leg)
    grid on
    print('-dpsc','-append','~/Documents/diurnalIsrael/results/winterOptimization.ps')
    saveas(gcf,'~/Documents/diurnalIsrael/results/winterOptimization.fig')
    % TODO - re run adfter the loop - and plot CC/CF for each case, winter and
    
end

if doSummer
    % summer optimization.
    
    %%
    % summer analysis
    close all
    disp('loading  summer results')
    load totalStationMat_summer_2to10.mat
    load optimizationWeights_summer_2to10
    optimizeFlag = 'summer'
    
    stationNumOptimizeVec = 2:10;
    
    if useOnlyLastResult
        stationNumOptimizeVec = stationNumOptimizeVec(end);
    end
    % running through optimization cases
    for stationNumOptimize=stationNumOptimizeVec
        simulationName = ['Israel optimized summer #',num2str(stationNumOptimize)];
        colors = jet(stationNumOptimize);
        
        % use optimization based stations
        clear weightedStations
        fig1 = figure;
        fig2 = figure;
        for i=1:stationNumOptimize
            if stationNumOptimize==1
                weightedStations = windEnergyNormalized(stations(1),:);
                weights{1} = 1;
            else
                weightedStations(i,:) = weights{stationNumOptimize}(i)*windEnergyNormalized(stations(i),:);
            end
            if find(tIsrael(stations(i),:)==0,1)
                loc = find(tIsrael(stations(i),:)==0,1);
                figure(fig1); hold on; plot(tIsrael(stations(i),1:loc-1),windEnergyNormalized(stations(i),1:loc-1),'.','color',colors(i,:))
                figure(fig2); hold on; plot(tIsrael(stations(i),1:loc-1),weightedStations(i,1:loc-1),'.','color',colors(i,:))
            else
                figure(fig1); hold on; plot(tIsrael(stations(i),:),windEnergyNormalized(stations(i),:),'.','color',colors(i,:))
                figure(fig2); hold on; plot(tIsrael(stations(i),:),weightedStations(i,:),'.','color',colors(i,:))
            end
        end
        [tVec,~,~,WindEnergyNormalizedAveraged, endVec, tfVec] = averageVectors(weightedStations,tIsrael(stations(1:stationNumOptimize),:),0);
        figure(fig1); hold on; plot(tVec,WindEnergyNormalizedAveraged,'k');
        datetick
        figure(fig1); legend({metaAll(stationVecOrig(stations(1:stationNumOptimize))).name,'averaged'})
        figure(fig2); hold on; plot(tVec,WindEnergyNormalizedAveraged,'k');
        datetick
        figure(fig2); legend({metaAll(stationVecOrig(stations(1:stationNumOptimize))).name,'averaged'})
        
        % re-run analysis for the weights{i} and stationVecOrig(stations(1:i))
        optimizationReRunAnalysisScript
        
        % plot map of stations and weights - as a circle of proportional size
        % to weights{i}(stationNum)
        figure(91);        
        disp('best staitions and weights --->')
        disp([stationVecOrig(stations(1:stationNumOptimize))' , stations(1:stationNumOptimize)', weights{stationNumOptimize}])
        plotStationLocations(stationVecOrig(stations(1:stationNumOptimize)),jet(83),weights{stationNumOptimize})

        print('-dpsc','-append',reportName)
        
        % calculate optimization parameter again
        if stationNumOptimize>1
            global input
            input.windEnergyNormalized = windEnergyNormalized(stations(1:stationNumOptimize),:);
            input.tIsrael = tIsrael(stations(1:stationNumOptimize),:);
            input.timeYearly = timeYearly;
            input.Prated = Prated;
            input.electricityNormalizedYearly = electricityNormalizedYearly;
            input.minNan = minNan;
            input.sortVector = sortVector;
            input.stationVecOrig = stations(1:stationNumOptimize);
            optimizationNormalized(stationNumOptimize) = -100*optimizeWeights(weights{stationNumOptimize})/CFIsrael(end-1); % [%] - reference year is 2011
        end
        
        % script for saving plots to report as .ps file
        reportDir = ['~/Documents/diurnalIsrael/results/reportSummer',num2str(stationNumOptimize),'/'];
        reportName = ['optimizationReport_summer_2to10_' num2str(stationNumOptimize) '.ps'];
        mkdir(reportDir)
        printOptimizationReport
        close all
        
    end
    
    % optimization convergence
    figure(1); hold on;
    plot(2:10,optimizationNormalized(2:10),'*-')
    title({'Minimum summer daily-peak capacity credit','as a function of stations used in optimization'})
    xlabel('Number of top ranking stations used')
    ylabel('Minimum summer daily-peak CC/CF [%]')
    ticks = get(gca,'YTick');
    for i=1:length(ticks)
        leg{i} = [num2str(ticks(i)),'%'];
    end
    set(gca,'YTickLabel',leg)
    axis tight
    grid on
    print('-dpsc','-append','~/Documents/diurnalIsrael/results/summerOptimization.ps')
    saveas(gcf,'~/Documents/diurnalIsrael/results/summerOptimization.fig')
    
    plotSimulatedPeakYear( maxDailyElectricity, dailyPeakCapacityFactor,yearlyCapacityFactor, stationName, year, stations, weights)
    
end