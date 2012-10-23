function plotPeak(stationName, Num, yearlyPeakCapacityFactorMinimum, yearlyPeakCapacityFactorMaximum, ...
                                    yearlyPeakCapacityFactorAverage, yearlyPeakCapacityFactorStd, ...
                                    monthlyPeakCapacityFactorMinimum, monthlyPeakCapacityFactorMaximum, ...
                                    monthlyPeakCapacityFactorAverage, monthlyPeakCapacityFactorStd, ...
                                    resultsDirectory, CF, years)
%plotPeak plots daily peak electricity demand statistics
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
if strcmp(stationName,'Israel wind')
    figure(890); hold on;
else if strcmp(stationName,'Israel solar')
        figure(891); hold on;
    else if strcmp(stationName,'Israel combined')
            figure(892); hold on;
        else
            figure(800+Num); hold on;
        end
    end
end
l = length(yearlyPeakCapacityFactorMinimum(Num,:));
if l>6 l=6; end
[h_ax,h_l] = plot2ylinscales(years(1:l),yearlyPeakCapacityFactorMinimum(Num,1:l),1/CF); % years may mismatch peakCapacityFactor
ylabel(h_ax(1),'CC')
ylabel(h_ax(2),'CC/CF')
set(h_l,'Marker','^'); set(h_l,'Color','m'); set(h_l,'LineStyle','-');set(h_l,'lineWidth',3);
plot(years(1:l), yearlyPeakCapacityFactorMaximum(Num,1:l),'c-v','lineWidth',3);
errorbar(years(1:l),yearlyPeakCapacityFactorAverage(Num,1:l),yearlyPeakCapacityFactorStd(Num,1:l)/2,'k-o','lineWidth',3)
title({[stationName 'daily peak capacity credit'],['Interannual capacity factor = ' num2str(100*CF,2) '%']});
a = legend('minimum','maximum','average and std');
set(get(a,'Title'),'String','Daily peak Capacity credit')
set(a,'Position',[0.3794 0.7057 0.2822 0.1335])
xlabel('year')
ylabel('Capacity credit')
set(h_ax(2),'Xtick',[])
set(h_ax(1),'Xtick',years)
axis tight
lin_fun = inline([num2str(1/CF) '*x']);
set(h_ax(2),'ylim',lin_fun(get(h_ax(1),'ylim')));
print([resultsDirectory 'YearlyPeakCC_', strrep(stationName,' ',''),'.png'],'-dpng')
saveas(gcf, [resultsDirectory 'YearlyPeakCC_', strrep(stationName,' ','')], 'fig')

% plotting monthly 
if strcmp(stationName,'Israel wind')
    figure(2490); hold on;
else if strcmp(stationName,'Israel solar')
        figure(2491); hold on;
    else if strcmp(stationName,'Israel combined')
            figure(2492); hold on;
        else
            figure(800+Num); hold on;
        end
    end
end
[h_ax,h_l] = plot2ylinscales(1:12,monthlyPeakCapacityFactorMinimum(Num,:),1/CF);
ylabel(h_ax(1),'CC')
ylabel(h_ax(2),'CC/CF')
set(h_l,'Marker','^'); set(h_l,'Color','m'); set(h_l,'LineStyle','-');set(h_l,'lineWidth',3);
plot(1:12,monthlyPeakCapacityFactorMaximum(Num,:),'c-v','lineWidth',3);
errorbar(1:12,monthlyPeakCapacityFactorAverage(Num,:),monthlyPeakCapacityFactorStd(Num,:)/2,'k-o','lineWidth',3)
title({[stationName ', interannual monthly-average, of daily peak capacity credit'],['Interannual capacity factor = ' num2str(100*CF,2) '%']});
a = legend('minimum','maximum','average and std');
grid on;
lin_fun = inline([num2str(1/CF) '*x']);
axis tight
axis([0 12 0 1])
set(h_ax(2),'ylim',lin_fun(get(h_ax(1),'ylim')));
set(h_ax(2),'Xtick',[])
set(get(a,'Title'),'String','Daily peak Capacity credit')
set(a,'Position',[0.1509    0.7109    0.2822    0.1335])
b = get(gcf,'Children');
set(b(3),'xTick',1:12)
set(b(3),'xTickLabel',months)
xlabel('month')
ylabel('Capacity credit')
print([resultsDirectory 'monthlyPeakCC_', strrep(stationName,' ',''),'.png'],'-dpng')
saveas(gcf, [resultsDirectory 'monthlyPeakCC_', strrep(stationName,' ','')], 'fig')
end

