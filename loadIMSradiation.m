
function [NIPDailyNorm_summer,NIPDailyNorm_winter,tDaily,ax,NIP,NIPDaily,d,m,y,h,mi] = loadIMSradiation(Num,fig,dataDirectory,p)
% SDED BOQER = 75
if ~(exist(dataDirectory))
    dataDirectory = '/home/hanan/Documents/measurements/';
end
% M = csvread('/VOLUMES/STUMPY32G/measurements/IMS-data/STATIONS DATA/SEDE BOQER/SEDE BOQER RADIATION/SEDE BOQER_RAD_2006.csv')
% open METADATA.csv (TODO - change to xls, from matlab) and look up station data
meta = loadMeta(Num);
pathname = [dataDirectory,'/IMS-data/STATIONS DATA/',meta.name,'/',meta.name,' RADIATION/'];
matFile = [pathname, 'Data_',num2str(Num),'.mat'];
if isempty(dir(matFile))
    % load data from station directory
    tic
    counter = 1;
    filenames = dir(pathname);
    for i=1:length(filenames)
        if or(~isempty(strfind(filenames(i).name,'csv')),~isempty(strfind(filenames(i).name,'CSV')))>0
            % read file
            DataFilename = [pathname filenames(i).name];
            disp(DataFilename)
            
            % read only data lines
            fid = fopen(DataFilename,'r');
            condition = 1;
            while condition
                Line = fgetl(fid);
                % check if end of file
                if Line == -1
                    break;
                end
                % jump over empty lines and header lines
                if and(~isempty(Line),~isempty(str2num(Line(1)))) 
                
                    %replace InVld with NaN
                    Line = strrep(Line,'InVld','NaN');
                    Line = strrep(Line,'Down','NaN');
                    Line = strrep(Line,'NoData','NaN');
                    Line = strrep(Line,'<Samp','NaN');
                    Line = strrep(Line,'Samp','NaN');
                    Line = strrep(Line,'Samp>','NaN'); 
                    % take out seconds statement if exists
                    loc = strfind(Line,':');
                    if length(loc)>1
                        Line(loc(2):(loc(2)+2))='';
                    end                 
                    dash = strfind(Line,'-');
                    % date is sometimes dd/mm/yy and sometimes dd-mm-yy           
                    if ~isempty(dash)
                        if dash(1)<4
                            Line(dash(1:2)) = '/'; %replace dd-mm-yy with dd/mm/yy and live - marks in the data intact
                        end
                    end
                    % read line
                    % data format
                    % 01/01/2006,16:50,-1,0,-3 - NIP,Grad,DiffR
                    [temp,num] = sscanf(Line,'%d/%d/%d,%d:%d,%g,%g,%g',[1,inf]);

                    d(counter) = temp(1); 
                    m(counter) = temp(2);
                    y(counter) = temp(3); 
                    h(counter) = temp(4); 
                    mi(counter) = temp(5);
                    NIP(counter) = temp(6);
                    Grad(counter) = temp(7);
                    DiffR(counter) = temp(8); 
                    disp(sprintf('%d/%d/%d %d:%d',y(counter),m(counter),d(counter),h(counter),mi(counter)));
                    counter = counter + 1;
                end
            end
            fclose(fid);
        end
    end    
    lengthTemp = i-1;
    if i==1
        disp(sprintf('problem with data'));
    end
    % saving workspace
    save(matFile)
    disp(['Reading ' DataFilename ' took ' num2str(toc) ' seconds'])
else
    load(matFile)
end

% diurnal plot for each month
dailyAvgTime = 1/6; % [hr]
tDaily = 0:dailyAvgTime:24;
M = length(tDaily)-1;
hTot = mi/60+h;     % hour time vector
monthString = {'January','February','March','April','May','June','July','August','September','October','November','December'};
col = jet(12);
for month=1:12
    for i=1:length(tDaily)-1
        loc = find(and(hTot>=tDaily(i),and(hTot<=tDaily(i+1),m==month)));
        NIPDaily(month,i) = nanmean(NIP(loc));
        GradDaily(month,i) = nanmean(Grad(loc));
        DiffRDaily(month,i) = nanmean(DiffR(loc));
        NIPDailyStd(month,i) = nanstd(NIP(loc));
        GradDailyStd(month,i) = nanstd(Grad(loc));
        DiffRDailyStd(month,i) = nanstd(DiffR(loc));
    end
end

NIPDailyNorm_summer = NIPDaily(7,:)/max(NIPDaily(7,:));
NIPDailyNorm_winter = NIPDaily(1,:)/max(NIPDaily(1,:));
zx = [];
if p
    figure(fig); hold on;
    ax = plot(tDaily(1:M)/24,NIPDailyNorm_summer,'r-.','lineWidth',4);
else
    ax=0;
end


%%
