function meta = loadMeta(N)

fid = fopen('METADATA.csv','r');

condition = 1; i = 1;
while condition
    str = fgetl(fid);
    if str == -1
        break;
    end   
    if ~isempty(str)
        metaStr = textscan(str,'%s','delimiter',';');
        %if nargin<1 disp(sprintf('%s\n',str)); end
        
        % num
        temp = strrep(metaStr{1}(1),'"','');
        num(i) = str2num(temp{1});
        % name (english)
        temp = strrep(metaStr{1}(3),'"','');
        name{i} = temp{1};
        % Israel coordinates
        IsrCoor1(i) = str2num(metaStr{1}{5});
        IsrCoor2(i) = str2num(metaStr{1}{6});
        % long-lat
        long(i) = str2num(metaStr{1}{7});
        lat(i) = str2num(metaStr{1}{8});
        % heights
        hAbsolute(i) = str2num(metaStr{1}{9});
        if ~isempty(str2num(metaStr{1}{11}))
            hGround(i) = str2num(metaStr{1}{11});
        else
            hGround(i) = -1; % designates no data
        end
        % is there even an anemometer installed? (if not - skip this data in plotIMS)
        if length(metaStr{1})<12 
            anemometer{i} = 'No';
        else
            anemometer{i} = metaStr{1}{12};
        end
    end
    i = i + 1; 
end

% find N in num
if (~exist('N'))
    for i=1:length(name)
        meta(i).num = num(i);
        meta(i).name = name{i};
        meta(i).IsrCoor(1) = IsrCoor1(i);
        meta(i).IsrCoor(2) = IsrCoor2(i);
        meta(i).long = long(i);
        meta(i).lat = lat(i);
        meta(i).h(1) = hAbsolute(i);
        meta(i).h(2) = hGround(i);
        meta(i).anemometer = anemometer{i};
    end
else
    loc = find(num==N);
    meta.num = num(loc);
    meta.name = name{loc};
    meta.IsrCoor(1) = IsrCoor1(loc);
    meta.IsrCoor(2) = IsrCoor2(loc);
    meta.long = long(loc);
    meta.lat = lat(loc);
    meta.h(1) = hAbsolute(loc);
    meta.h(2) = hGround(loc);
    meta.anemometer = anemometer{loc};
end
fclose(fid);