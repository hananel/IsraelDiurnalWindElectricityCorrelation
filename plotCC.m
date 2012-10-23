function plotCC(percentVec, sp1,sp2,Num,year,yearLegend,CFper, station,...
                CFperIsraelS, CFperIsraelSW)
%plotCC plots capacity credit as function of precentage of peak electricity
%use
%   plotCC plots capacity credit for a multi year analysis, representing each year seperately in a
%   subplot
                             
                       

if not(station)
    figure(10015);
    subplot(sp1,sp2,yearLegend);hold on;
    percentCounter = 0;
    for i=percentVec
        percentCounter = percentCounter + 1;
        CF_temp(percentCounter) = CFper(percentCounter);
    end    
else
    figure(100+Num);  
    subplot(sp1,sp2,yearLegend);hold on;
    percentCounter = 0;
    for i=percentVec
        percentCounter = percentCounter + 1;
        CF_temp(percentCounter) = CFper(1,1,percentCounter);
    end
end
CF = CF_temp(end); 
[h_ax,h_l] = plot2ylinscales(percentVec,CF_temp,1/CF);
ylabel(h_ax(1),'CC')
ylabel(h_ax(2),'CC/CF')
set(h_ax(2),'Xtick',[])
set(h_l,'Color','b');
if not(station) % for checking all of Israel or a few stations vs. solar as well as wind
    plot(percentVec,CFperIsraelS,'r'); 
    plot(percentVec,CFperIsraelSW,'g');
    maxCF = 1;
else
    maxCF = nanmax(CF_temp);  
end
title([num2str(year) ' CF=' num2str(100*CF,2) '%']);
set(gca,'xtick',[0,30, 50, 100]);
minCF = 0;
axis tight
lin_fun = inline([num2str(1/CF) '*x']);
set(h_ax(2),'ylim',lin_fun(get(h_ax(1),'ylim')));
% plotting 30% and 50% lines          
plot([30 30],[0 1],'b')
plot([50 50],[0 1],'r')

end

