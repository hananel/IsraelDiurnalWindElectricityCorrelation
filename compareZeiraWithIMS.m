function compareZeiraWithIMS()
% compares diurnal pattern of zeira measurements at heights 30,40,50,60 with pichman 4 km away at 10 meter

% lat/lon gof katif station (zeira data)
lat = 33.038691; long = 35.833517;
stationNum = 1;


global dataDirectory
dataDirectory = '/home/hanan/Documents/measurements/'
global resultsDirectory
resultsDirectory = '/home/hanan/Documents/measurements/';
debug_on_warning(0);
debug_on_error(0);
more off;
close all;

% wind rose difinition
roseDir = [0:10:350]';

% 1
% 
% finding closest 3 stations and plotting their location and distance from the point

% loading meta
meta = loadMeta;

% plotting all stations
figure(100);

% zoom out plot - all of Israel
subplot(2,2,[1,3]);
plotIsraelBorders
axis tight; hold on;
g = plot(long,lat,'or','markersize',6);
set(g,'markerfacecolor','r')
ax = get(gcf); 
set(ax.children,'position',[0.05000   0.11000   0.43466   0.81500]);
set(ax.children,'xlim',[33.5   35.892]);
for i=1:length(meta)
    % checking for anemometer
    yesAnemometer = not(or(strcmp(meta(i).anemometer,'No'),meta(i).h(2)==-1));
    if yesAnemometer
        text(meta(i).long,meta(i).lat,num2str(meta(i).num),'color','g');
        % building distance vector
        distance(i) = haversine(meta(i).lat,meta(i).long,lat,long);
    else
        text(meta(i).long,meta(i).lat,num2str(meta(i).num),'color','b');
        distance(i) = 999;
    end

end    

% sorting to find closest stations
distanceOrig = distance;
[distance,stations] = sort(distance);

for i=1:stationNum
    text(meta(stations(i)).long,meta(stations(i)).lat,num2str(meta(stations(i)).num),'color','r','fontweight','bold')
    plot([long meta(stations(i)).long],[lat meta(stations(i)).lat],'r');
end 
xlabel('longitude [deg]'); ylabel('latitude [deg]'); title('IMS stations');

% zoom in plot
zoomDistance = 0.3; % [deg]
ax = subplot(2,2,4);
plotIsraelBorders
axis tight; hold on;
g = plot(long,lat,'or','markersize',6);
set(g,'markerfacecolor','r')
for i=1:length(meta)
    % checking if the station is within the zoom area
    yesZoom = distanceOrig(i)<(0.75*haversine(lat,long,lat+zoomDistance,long+zoomDistance));
    % checking for anemometer
    yesAnemometer = not(or(strcmp(meta(i).anemometer,'No'),meta(i).h(2)==-1));
    if yesZoom
        if yesAnemometer
            text(meta(i).long,meta(i).lat,num2str(meta(i).num),'color','g');
        else
            text(meta(i).long,meta(i).lat,num2str(meta(i).num),'color','b');
        end
    end
end 

% finding the first stationNum stations
for i=1:stationNum
    text(meta(stations(i)).long,meta(stations(i)).lat,num2str(meta(stations(i)).num),'color','r','fontweight','bold')
    plot([long meta(stations(i)).long],[lat meta(stations(i)).lat],'r');
end 
xlabel('longitude [deg]'); ylabel('latitude [deg]'); title('Zoom in on close stations');
axis([long-zoomDistance long+zoomDistance lat-zoomDistance lat+zoomDistance])

% report of location of closest station
subplot(2,2,2);
axis([0 1 0 1],'off');
text(-0.1,0.9,['The location of the building is ',num2str(long,4),'/',num2str(lat,4)]);
text(-0.1,0.8,['The closest ',num2str(stationNum),' meteorological stations with ']);
text(-0.1,0.7,'long term wind data are:');
for i=1:stationNum
    text(-0.1,0.6-i*0.2,[num2str(meta(stations(i)).num),': ',meta(stations(i)).name,' ', num2str(distance(i),2),' km away']);
    text(-0.1,0.6-i*0.2-0.1,['with anemometer at ', num2str(meta(stations(i)).h(2)), ' meter']);
end

% print report
reportDirectory = ['lat_',num2str(lat,4),'_long_',num2str(long,4)];
mkdir(resultsDirectory,reportDirectory);
print([resultsDirectory,'/',reportDirectory,'/IMS_stations_near_long_', num2str(long,4),'_lat_',num2str(lat,4),'.png']);

% 2
%
% loading station data and submitting report for each
% TODO - at the moment this is for the average day of each month. should be representative enough, but I should compare to IMS publications.
disp(['report for ',num2str(stationNum),' stations'])
metaAll = meta;
Num = stations(1);
legendText = {};
col = jet(stationNum);
% load data
pathname = [dataDirectory,'/IMS-data/STATIONS DATA/',metaAll(Num).name,'/',metaAll(Num).name,'/'];
matFile = [pathname, 'Data_',num2str(Num),'.mat'];
disp('loading matFile')
% clearing UMonthly
if exist('Umonthly')
    clear Umonthly;
end
load(matFile);
M = length(tDaily)-1;
% make station directory
stationDirectory = strrep(meta.name,' ','');
mkdir([resultsDirectory,'/',reportDirectory],stationDirectory); 
legendText = {legendText{:},meta.name}

% calculating average yearly data - by a running average filter
% first connecting the month time line
%UYear = []; 
for i=1:12
    %UYear = [UYear,UMonthly(i,:)];
    um(i) = nanmean(U(find(m==i)));
    ustd(i)=nanstd(U(find(m==i)));
end

% calculating diurnal for only 2011
if ~exist('UDaily2011')
    for month=1:12
        for i=1:length(tDaily)-1
            % 2011
            loc = find(and(hTot>=tDaily(i),hTot<=tDaily(i+1),m==month,y==2011));
            UDaily2011(month,i) = nanmean(U(loc));
            UStdDaily2011(month,i) = nanstd(U(loc));
        end
    end
end

M = length(tDaily)-1;
totalHours = 0;
tDailyS = tDaily; % name substituting - not to run over with the "load" action
MS = M;

% loading Zeira data
ZeiraMatFile = '/home/hanan/Dropbox/MyPHDArticles/InProgress/IsraelDiurnalWind/zeira/Data_zeira.mat';
load(ZeiraMatFile);

sprintf('average wind speed is %g\n',nanmean(U))
for month=1:12
    
    % non normalized plots
    figure(month+10); 
    subplot(211)
    ax = plot(tDailyS(1:MS)/24,UDaily(month,1:MS),'*');
    hold on;
    set(ax,'color',col(month,:))
    set(gca,"xtick",[0,0.25,0.5,0.75,1])
    title(sprintf(['Wind speed [m/s] from ',meta.name,' compared with Katif for ', monthString{month}]));
    ax = plot(tDaily(1:M)/24,UDaily30(month,1:M),'o');
    set(ax,'color',col(month,:))
    ax = plot(tDaily(1:M)/24,UDaily40(month,1:M),'^');
    set(ax,'color',col(month,:))
    ax = plot(tDaily(1:M)/24,UDaily50(month,1:M),'v');
    set(ax,'color',col(month,:)) 
    ax = plot(tDaily(1:M)/24,UDaily60(month,1:M),'+');
    set(ax,'color',col(month,:))
    xlabel('Hour'); ylabel('U [m/s]'); 
    axis([0,1,0,nanmax(nanmax(UDaily)+nanmax(UStdDaily))])
    datetick('x',15,'keeplimits','keepticks');
    legend({[meta.name,' 10 m'],'Katif 30 m', 'Katif 40 m','Katif 50 m','Katif 60 m'})
    
    % normalized energy plots
    subplot(212) 
    ax = plot(tDailyS(1:MS)/24,UDaily(month,1:MS).^3/max(UDaily(month,1:MS).^3),'*');
    hold on;
    set(ax,'color',col(month,:))
    set(gca,"xtick",[0,0.25,0.5,0.75,1])
    title(sprintf(['Normalized wind energy from ',meta.name,' compared with Katif for ', monthString{month}]));
    ax = plot(tDaily(1:M)/24,UDaily30(month,1:M).^3/max(UDaily30(month,1:M).^3),'o');
    set(ax,'color',col(month,:))
    ax = plot(tDaily(1:M)/24,UDaily40(month,1:M).^3/max(UDaily40(month,1:M).^3),'^');
    set(ax,'color',col(month,:))
    ax = plot(tDaily(1:M)/24,UDaily50(month,1:M).^3/max(UDaily50(month,1:M).^3),'v');
    set(ax,'color',col(month,:)) 
    ax = plot(tDaily(1:M)/24,UDaily60(month,1:M).^3/max(UDaily60(month,1:M).^3),'+');
    set(ax,'color',col(month,:))
    xlabel('Hour'); ylabel('Normalized energy'); 
    axis([0,1,0,1])
    datetick('x',15,'keeplimits','keepticks');
    legend({[meta.name,' 10 m'],'Katif 30 m', 'Katif 40 m','Katif 50 m','Katif 60 m'},'location','south')
    
    % plotting the maximum of each curve
    [max10_2011,imax10_2011(month)] = max(UDaily2011(month,1:MS));
    [max10,imax10(month)] = max(UDaily(month,1:MS));
    [max30,imax30(month)] = max(UDaily30(month,1:M));
    [max40,imax40(month)] = max(UDaily40(month,1:M));
    [max50,imax50(month)] = max(UDaily50(month,1:M));
    [max60,imax60(month)] = max(UDaily60(month,1:M));
    plot([tDailyS(imax10(month)) tDailyS(imax10(month))]/24,[0 1],'r')
    plot([tDaily(imax30(month)) tDaily(imax30(month))]/24,[0 1],'b')
    plot([tDaily(imax40(month)) tDaily(imax40(month))]/24,[0 1],'b')
    plot([tDaily(imax50(month)) tDaily(imax50(month))]/24,[0 1],'b')
    plot([tDaily(imax60(month)) tDaily(imax60(month))]/24,[0 1],'b')
    
    % print
    print([resultsDirectory,'/Diurnal_energy_month_', monthString{month},'.png']);
    
end


% plotting maximum time change
months = {'January','February','March','April','May','June','July','August','September','October','November','December'}; 
measVec = [10 30 40 50 60];
figure(555);
badM = [9 11 12 2 3];

for month=1:12
    if sum(month == badM)
        ax = plot([tDailyS(imax10(month)) tDaily([imax30(month),imax40(month),imax50(month),imax60(month)])]/24,measVec,'-o');
    else
        ax = plot([tDailyS(imax10(month)) tDaily([imax30(month),imax40(month),imax50(month),imax60(month)])]/24,measVec,'-o','lineWidth',3);
    end
    set(ax,'color',col(month,:))
    hold on;
end
legend(months)
title('maximum time of wind energy vs. measurement height. Picman (10 m) and Katif (30-60 m)')
xlabel('Hour'); ylabel('height [m]');
set(gca,"xtick",[0,0.25,0.5,0.75,1])
datetick('x',15,'keeplimits','keepticks');

% print
print([resultsDirectory,'/Diurnal_energy_maxHourWithHeight_', monthString{month},'.png']);
    
% plotting maximum time change - with x as month vec
monthNum = [datenum(0,1,15),datenum(0,2,15),datenum(0,3,15), ...
                datenum(0,4,15),datenum(0,5,15),datenum(0,6,15), ...
                datenum(0,7,15),datenum(0,8,15),datenum(0,9,15), ...
                datenum(0,10,15),datenum(0,11,15),datenum(0,12,15)];
figure(556);
ax = plot(monthNum,tDailyS(imax10)/24,'k-*'); hold on;
ax = plot(monthNum,tDaily(imax30)/24,'r^');
ax = plot(monthNum,tDaily(imax40)/24,'bv');
ax = plot(monthNum,tDaily(imax50)/24,'g+');
ax = plot(monthNum,tDaily(imax60)/24,'mo');
hold on;
legend({'10','30','40','50','60'},'location','northwest')
set(gca,'xtick',monthNum)
axis([monthNum(1)-monthNum(1)/2 monthNum(end)+monthNum(1)/2 0 1])
datetick('x',4,'keeplimits','keepticks');
xlabel('months'); ylabel('Hour');
title('maximum time of wind energy vs. measurement height. Picman (10 m) and Katif (30-60 m)')
set(gca,"ytick",[0,0.25,0.5,0.75,1])
datetick('y',15,'keeplimits','keepticks');

% print
print([resultsDirectory,'/Diurnal_energy_maxHourWithHeight_perMonth_', monthString{month},'.png']);
    
figure(557)
ax = plot(monthNum,tDailyS(imax10_2011)/24,'k-*'); hold on;
ax = plot(monthNum,tDaily(imax30)/24,'r^');
ax = plot(monthNum,tDaily(imax40)/24,'bv');
ax = plot(monthNum,tDaily(imax50)/24,'g+');
ax = plot(monthNum,tDaily(imax60)/24,'mo');
hold on;
legend({'10-2011','30','40','50','60'},'location','northwest')
set(gca,'xtick',monthNum)
axis([monthNum(1)-monthNum(1)/2 monthNum(end)+monthNum(1)/2 0 1])
datetick('x',4,'keeplimits','keepticks');
xlabel('months'); ylabel('Hour');
title('maximum time of wind energy vs. measurement height. Picman-2011 (10 m) and Katif (30-60 m)')
set(gca,"ytick",[0,0.25,0.5,0.75,1])
datetick('y',15,'keeplimits','keepticks');

% print
print([resultsDirectory,'/Diurnal_energy_maxHourWithHeight_perMonth_2011_', monthString{month},'.png']);
  
function d = haversine(lat1,lon1,lat2,lon2)
% a = sin²(Δlat/2) + cos(lat1).cos(lat2).sin²(Δlong/2)
% c = 2.atan2(√a, √(1−a))
% d = R.c
R = 6371; % km
dLat = (lat2-lat1)*pi/180;
dLon = (lon2-lon1)*pi/180;
lat1 = lat1*pi/180;
lat2 = lat2*pi/180;
a = (sin(dLat/2))^2 + (sin(dLon/2))^2 * cos(lat1) * cos(lat2); 
c = 2 * atan2(sqrt(a), sqrt(1-a)); 
d = R * c;

function makeReport(s)
% accepts structure s and creates pdf report and png graphs
