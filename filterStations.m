function [stationVec, colOld]  = filterStations( stationVec, Ulimit)
%filterStations filters stationVec according to minimum Ulimit ([m/s]) and
% low wind energy locations
stationVecOld = stationVec;
UmeanVec = csvread('stationUmean.csv');
if stationVec>0 % user wants a specific set of stations
	stationVec = stationVecOld;
	stationVec = ismember(stationVec,find(round(UmeanVec*10)/10>=Ulimit)).*stationVec;
	stationVec(stationVec==0)=[];
	% throwing away zefat har kenan - the peak of july is very different from the rest, doesn't make sense
	stationVec(stationVec==9)=[];
	% throwing away coast stations
    stationVec(stationVec==5)=[];
    stationVec(stationVec==8)=[];
    stationVec(stationVec==19)=[];
    stationVec(stationVec==21)=[];
    stationVec(stationVec==29)=[];
	stationVec(stationVec==31)=[];
    stationVec(stationVec==44)=[];
	stationVec(stationVec==48)=[];
	stationVec(stationVec==62)=[];
	stationVec(stationVec==35)=[];
	% throwing away bikaa stations
	stationVec(stationVec==74)=[];
	stationVec(stationVec==51)=[];
	stationVec(stationVec==64)=[];
	stationVec(stationVec==66)=[];
else        % filtering all IMS stations 
	% finding the stations that are above Ulimit
	stationVec = find(round(UmeanVec*10)/10>=Ulimit);
	% throwing away zefat har kenan - the peak of july is very different from the rest, doesn't make sense
	stationVec(stationVec==9)=[];
	% throwing away coast stations
    stationVec(stationVec==5)=[];
    stationVec(stationVec==8)=[];
    stationVec(stationVec==19)=[];
    stationVec(stationVec==21)=[];
    stationVec(stationVec==29)=[];
    stationVec(stationVec==31)=[];
	stationVec(stationVec==44)=[];
	stationVec(stationVec==48)=[];
	stationVec(stationVec==62)=[];
	stationVec(stationVec==35)=[];
	% throwing away bikaa stations
	stationVec(stationVec==74)=[];
	stationVec(stationVec==51)=[];
	stationVec(stationVec==64)=[];
	stationVec(stationVec==66)=[];
	% adding Eilat station
	stationVec(end+1) = 83;
end
disp(sprintf('\nTotal number of stations used is %d',length(stationVec)));
colOld = jet(length(stationVec)); 
figure(1); 	
plotStationLocations(stationVec,colOld);
end

