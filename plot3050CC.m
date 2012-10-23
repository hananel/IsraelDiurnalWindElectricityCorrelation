function plot3050CC(Num, stationName,resultsDirectory, windDataYearVec,CFtemp,CF50,CF30,CF1,Cor,yearVec)
%plot3050CC Plots the 30%, 50% and correlation between a stations
%yearly wind energy production and average electricity yearly consumption
%trend 

CF = nanmean(CFtemp);
figure(103+Num); hold on;
title({'Yearly 1/30/50 capacity credits vs. correlation',['Average interannual capacity factor = ' num2str(100*CF,2) '%']})
set(gcf,'name',stationName)
plot(windDataYearVec,CF50,'r-v',windDataYearVec,CF30,'b-+',windDataYearVec,CF1,'c-+', windDataYearVec,Cor,'g-o');
[h_ax,h_l] = plot2ylinscales(windDataYearVec,CFtemp,1/CF);
ylabel(h_ax(1),'CC, CF, corrlation')
ylabel(h_ax(2),'CC/CF')
set(h_l,'Marker','^'); set(h_l,'Color','k');
axis tight
lin_fun = inline([num2str(1/CF) '*x']);
set(h_ax(2),'ylim',lin_fun(get(h_ax(1),'ylim')));
legend('CC_{50}','CC_{30}','CC_{1}','correleation','CF','Location','Best');
xlabel('year'); title([stationName ', interannual capacity factor = ' num2str(100*CF,2) '%'])
set(gca,'xtick',yearVec);
set(h_ax(2),'Xtick',[])
grid on;
print([resultsDirectory 'CC_', strrep(stationName,' ',''),'.png'],'-dpng');           
saveas(gcf, [resultsDirectory 'CC_', strrep(stationName,' ','')], 'fig')
end

