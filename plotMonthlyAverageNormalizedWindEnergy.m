function plotMonthlyAverageNormalizedWindEnergy( Num, stationName, electricityNormalizedDaily, UDaily, months, resultsDirectory )
%plotMonthlyAverageNormalizedWindEnergy plots normalized wind energy (not
%turbine

figure(203+Num); hold on;
set(gcf,'name',stationName)
for month=1:12
    subplot(4,3,month); hold on;
    plot(0:23,electricityNormalizedDaily(month,:),'k','LineWidth',3);
    EDailyNorm(month,:) = UDaily(month,:).^3/max(UDaily(month,:).^3);
    plot((0:143)/6,EDailyNorm(month,:),'r--','LineWidth',3)
    title(months(month))
    axis([0 24 0 1])
    if sum(month==[2,3,5,6,8,9])
        axis off
    else if sum(month==[1,4,7])
             set(gca,'xtick',[-100])
             set(gca,'ytick',[0,1]);
         else if sum(month==[11,12])
                set(gca,'xtick',[0,12,24]);
                set(gca,'ytick',[-100])
             else
                set(gca,'xtick',[0,12,24]);
                set(gca,'ytick',[0,1]);
             end
         end
    end
end
set(gcf,'Color',[1 1 1])
subplot(4,3,11); xlabel(stationName);
print([resultsDirectory 'monthlyWindEnergyCorrelation_', strrep(stationName,' ',''),'.png'],'-dpng');
saveas(gcf, [resultsDirectory 'monthlyWindEnergyCorrelation_', strrep(stationName,' ','')], 'fig')
end

