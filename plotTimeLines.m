function plotTimeLines(D,tD,endVec,tfVec,col,M,tM, stationNum)

global debug
%maximum time range
ts = min(tD(:,1));
tf = max(tfVec);
if nargin>4
    subPlotNum = min(size(tD)+1);
    if debug
        fig = figure(2000);
        for i=1:subPlotNum-1
            subplot(subPlotNum,1,i); hold on;
            plot(tD(i,1:endVec(i)),D(i,1:endVec(i)),'.','color',col(i,:))
            axis([ts tf 0 max(max(D))])
            axis off
        end
        subplot(subPlotNum,1,subPlotNum);
        plot(tM,M,'.k')
        axis([ts tf 0 max(max(D))])
        datetick('x',10,'keeplimits','keepticks');
        box off
        set(gca,'ycolor','w')
        set(gca,'ytick',[])
        % set(gca,
        print('TimeLinePlots.png','-dpng')
    end
	fig = figure(3000); hold on;
    title({'IMS stations used' 'Available wind speed data'})
	for i=1:(subPlotNum-1)
        %TODO
		%plot(tD(i,1:endVec(i)),i*not(isnan(D(i,1:endVec(i)))),'.','lineWidth',2,'color',col(i,:))
		plot(tD(i,1:endVec(i)),(subPlotNum-i)*(ones(1,endVec(i))),'.','lineWidth',2,'color',col(i,:))
        plot([tD(i,1) tD(i,1)],[subPlotNum-i-0.25 subPlotNum-i+0.25],'lineWidth',2,'color',col(i,:))
		plot([tD(i,endVec(i)) tD(i,endVec(i))],[subPlotNum-i-0.25 subPlotNum-i+0.25],'lineWidth',2,'color',col(i,:))
	end
	plot(tM,subPlotNum*ones(1,length(tM)),'.k')
	plot([tM(1) tM(1)],[1-0.25 subPlotNum+0.25],'k')
	plot([tM(end) tM(end)],[1-0.25 subPlotNum+0.25],'k')
	axis([ts tf 0 subPlotNum+1])
	datetick('x',10,'keeplimits');
    yy = get(gca,'YTick');
    set(gca,'YTick',1:size(tD,1)+1)
	set(gca,'YTickLabel',{stationNum(end:-1:1),'Total'})
	box off
	print('TimeLinePlots_simple.png','-dpng')
	% TODO - another way to check the results visualy?
else
	fig = figure;
	subPlotNum = min(size(tD));
	for i=1:subPlotNum
		subplot(subPlotNum,1,i);
		maxt = find(tD(i,:)==0,1);
		if isempty(maxt) maxt = length(tD(i,:)); end
		plot(tD(i,1:maxt-1),D(i,1:maxt-1),'.')
		datetick('x',10,'keeplimits','keepticks');
	end
end
