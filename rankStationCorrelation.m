function [bestStationSetNum, bestStationSetVec, sortedStationVec, bestStationSetWeight, finalRank ,a ,b] = rankStationCorrelation(CorDaily, Cor, stationVecOrig,season, metaAll)
%rankStationsCorrelation calculate station ranking based on correlation so that I can try to
% cross correlt the stations themselves!
% basis - look for stations with most days that are non correlated (look at
% correlation of the electricity-correlation matrix

global debug
global report

% find station-correlation between daily electricity-correlation
disp('calculating correlation between all stations')
colors = jet(365);

% throw away nan from corDaily
CorDaily(isnan(CorDaily)) = 0;

yearStart = 5; %2010
yearEnd   = 6; %2011
for year=yearStart:yearEnd
    temp = CorDaily(stationVecOrig,year,:);
    temp = reshape(temp,[size(temp,1),365]);
    if debug
        figure(99); 
        bar(temp')
        xlabel('day')
        ylabel('daily correlation')
        title(year)
        saveas(gcf,['rankingDailyCorrelation_26Stations_',num2str(year),'.fig'])
    end
    disp(year)
    
    switch season
        case 'all'
            if year==1
                disp('according to all year daily corellations'); end
            % looking for best combination of a given number of stations -
            % sorting by daily correlation for the entire year
            [tempSorted,tempSortVec] = sort(temp,'descend');
            % rank the best choices
            for i=1:length(stationVecOrig)
                % take the sorted stations - and save the number one choice, and
                % any subsequent choice (TODO - the 2nd choice should be the first
                % + the new added stations, and so on.
                sortedStations = sort(unique(tempSortVec(i,:)));
                % find number of occurance of each station in the "top" list
                clear occuranceOfStation
                for j=1:length(sortedStations)
                    occuranceOfStation(j) = sum(tempSortVec(i,:)==sortedStations(j));
                end
                % now sort this - putting up front the stations with most days of
                % best correlation
                [occurance{year,i},bestStationSetVec{year,i}] = sort(occuranceOfStation,'descend');
                bestStationSetNum{year,i} = stationVecOrig(sortedStations(bestStationSetVec{year,i}));
                sortedStationVec{year,i} = sortedStations;
                
                %now - need to weight the capacity of the stations.
                % options:
                % 1. by the ratio of the aquared correlation?
                % 2. by the best correlated days count (same as used for ranking the
                % stations themselves)
                % TODO - now leaving at one/length(stationVecOrig)
                bestStationSetWeight{year,i} = ones(1,length(sortedStations)) / length(sortedStations);
            end
        case 'winter'
            if year==1
                disp('according to winter daily corellations'); end
            %looking for best daily combination for winter
            % winter days - 1-70, 335-365
            loc = [1:70,335:365];
            [tempSorted,tempSortVec] = sort(temp(:,loc),'descend');
            % rank the best choices
           
            for i=1:length(stationVecOrig)
                % take the sorted stations - and save the number one choice, and
                % any subsequent choice (TODO - the 2nd choice should be the first
                % + the new added stations, and so on.
                sortedStations = sort(unique(tempSortVec(i,:)));
                % find number of occurance of each station in the "top" list 
                clear occuranceOfStation
                for j=1:length(sortedStations)
                    occuranceOfStation(j) = sum(tempSortVec(i,:)==sortedStations(j));
                end
                % now sort this - putting up front the stations with most days of
                % best correlation
                [occurance{year,i},bestStationSetVec{year,i}] = sort(occuranceOfStation,'descend');
                bestStationSetNum{year,i} = stationVecOrig(sortedStations(bestStationSetVec{year,i}));
                sortedStationVec{year,i} = sortedStations;
                
                %now - need to weight the capacity of the stations.
                % options:
                % 1. by the ratio of the aquared correlation?
                % 2. by the best correlated days count (same as used for ranking the
                % stations themselves)
                % TODO - now leaving at one/length(stationVecOrig)
                bestStationSetWeight{year,i} = ones(1,length(sortedStations)) / length(sortedStations);
            end
        case 'summer'
            if year==1
                disp('according to winter daily corellations'); end
            %looking for best daily combination for summer
            % summer days - 150-340
            loc = 170:250;
            [tempSorted,tempSortVec] = sort(temp(:,loc),'descend');
            % rank the best choices
            for i=1:length(stationVecOrig)
                % take the sorted stations - and save the number one choice, and
                % any subsequent choice (TODO - the 2nd choice should be the first
                % + the new added stations, and so on.
                sortedStations = sort(unique(tempSortVec(i,:)));
                % find number of occurance of each station in the "top" list
                clear occuranceOfStation
                for j=1:length(sortedStations)
                    occuranceOfStation(j) = sum(tempSortVec(i,:)==sortedStations(j));
                end
                % now sort this - putting up front the stations with most days of
                % best correlation
                [occurance{year,i},bestStationSetVec{year,i}] = sort(occuranceOfStation,'descend');
                bestStationSetNum{year,i} = stationVecOrig(sortedStations(bestStationSetVec{year,i}));
                sortedStationVec{year,i} = sortedStations;
                
                %now - need to weight the capacity of the stations.
                % options:
                % 1. by the ratio of the aquared correlation?
                % 2. by the best correlated days count (same as used for ranking the
                % stations themselves)
                % TODO - now leaving at one/length(stationVecOrig)
                bestStationSetWeight{year,i} = ones(1,length(sortedStations)) / length(sortedStations);
            end
    end
end

% return interannual best match (only number 1 ranks at the moment), by weighting them according to each year's
% rank
finalRank = zeros(1,length(stationVecOrig));
for year=yearStart:yearEnd
    for i=1:length(bestStationSetNum{year,1})
        if ~isempty(find(bestStationSetNum{year,i}==stationVecOrig(i)))
            finalRank(i) = finalRank(i) + find(bestStationSetNum{year,i}==stationVecOrig(i));
        end
    end
end

if report
    for i =1:15
        b(i,1:length(occurance{yearStart,i})) = occurance{yearStart,i};
    end
    for i =1:15
        a(i,1:length(bestStationSetNum{yearStart,i})) = bestStationSetNum{yearStart,i};
    end
    a(a==0) = NaN;
    b(b==0) = NaN;
    
    figure;
    barh(a',b')
    legend('1st','','','','','','','','','','','','','','15th')
    stations = unique(unique(a));
    stations(isnan(stations)) = [];
    set(gca,'YTick',stations)
    %hText = xticklabel_rotate(get(gca,'XTick'),45,{metaAll(stations).name});
    set(gca,'YTickLabel',{metaAll(stations).name})
    set(gca,'YDir','reverse')
    title(['stations best correlated to electricity load during ',season ' seasons'])
    xlabel('days of best correlation to electricity load')
else
    a = 0; b = 0;
end
end

