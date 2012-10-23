function [electricityNormalizedYearly,tYearly,electricityNormalizedDaily] = averageIECyear(dataPath,plotme,years)
% averageIECyear averages the IEC forcast data to an average year and saves
% csv files for [electricityNormalizedYearly,tYearly,electricityNormalizedDaily]

% TODO - save csv files of
% [electricityNormalizedYearly,tYearly,electricityNormalizedDaily] after
% you change their names!

if sum(dataPath==0)==0
    load([dataPath 'load_base_case.mat']);
    t = double(t);
    w = double(w);
    timeLine=datenum(t(:,1),t(:,2),t(:,3),t(:,4),zeros(length(t(:,1)),1),zeros(length(t(:,1)),1));
    
    % normalizing each year of data
    counter = 1; wNorm = double(w);
    wYearlyMax = ones(1,length(unique(t(:,1))));
    iVec = unique(t(:,1));
    yearVec = iVec;
    
    for counter=1:length(iVec)
        loc = find(t(:,1)==iVec(counter));
        wYearlyMax(counter) = double(max(w(loc)));
        wNorm(loc) = double(w(loc))/wYearlyMax(counter);
    end

    % breaking into diurnal average per month
    months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
    col = jet(12);
    if nargin<3
        years = yearVec;
    end
    
    for month=1:12
        for h=0:23
            loc = find(and(t(:,4)>=h,and(t(:,4)<=h+1,and(t(:,2)==month,ismember(t(:,1),years)))));
            Edaily(month,h+1) = double(nanmean(w(loc)));
            electricityNormalizedDaily(month,h+1) = double(nanmean(wNorm(loc)));
        end
    end

    % averaging - to get a representative year
    counter=1;
    for month=1:12
        for day=1:max(t(find(t(:,2)==month),3))
            for h=0:23
                loc = find(and(and(t(:,2)==month,ismember(t(:,1),years)),and(t(:,3)==day,t(:,4)==h)));
                electricityNormalizedYearly(counter) = nanmean(wNorm(loc));
                counter = counter + 1;
            end
        end
    end
    % final normalizing - so that yearly peak is 1
    electricityNormalizedYearly = electricityNormalizedYearly / max(electricityNormalizedYearly);
    
    % saving
    tYearly = 0:length(electricityNormalizedYearly)-1;
    csvwrite('averageIECyearly.csv',[tYearly' electricityNormalizedYearly'])
    csvwrite('averageIECdaily.csv',electricityNormalizedDaily);
else
    electricityNormalizedDaily = csvread('averageIECdaily.csv');
    temp = csvread('averageIECyearly.csv');
    tYearly = temp(:,1);
    electricityNormalizedYearly = temp(:,2);
end

if nargin>1
    figure(11); subplot(211)
    plot(timeLine,w,'.'); 
    datetick('x','keeplimits','keepticks');
    ylabel('Mw')
     figure(11); subplot(212)
    plot(wNorm,'.'); 
    %datetick('x','keeplimits','keepticks');
    ylabel('Yearly normalized Mw')
    
    figure();
    for month=1:12    
        plot([0:23]/24,electricityNormalizedDaily(month,:),'color',col(month,:));
        hold on;
    end
    xlabel('Hour'); ylabel('Normalized consumed electricity'); 
    axis([0,1,0,1])
    set(gca,'xtick',[0,0.25,0.5,0.75,1])
    datetick('x',15,'keeplimits','keepticks');
    legend(months,'location','south')  
    if length(years)==1
        title({'averaged and normalized representative electricity diurnal patterns - for forecast data',['year: ' num2str(years(1))]})
    else
        title({'averaged and normalized representative electricity diurnal patterns - for forecast data',['years: ' num2str(years(1)) , ' to ', num2str(years(end))]})
    end
    figName = ['IEC_diurnal4_',num2str(years)];
    print([figName '.eps'],'-depsc')
    saveas(gcf,figName,'fig');
    
    figure;
    plot(tYearly/365/24*12,electricityNormalizedYearly); 
    if length(years)==1
        title({'averaged and normalized representative electricity year - for forecast data',['year: ' num2str(years(1))]})
    else
        title({'averaged and normalized representative electricity year - for forecast data',['years: ' num2str(years(1)) , ' to ', num2str(years(end))]})
    end
    figName = ['IEC4_',num2str(years)];
    print([figName '.eps'],'-depsc')
    saveas(gcf,figName,'fig');
end


