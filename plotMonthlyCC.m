function plotMonthlyCC(stationName, CFtempMonthly, CF50Monthly, CF30Monthly, CF1Monthly, CF)
%plotMonthlyCF plot monthly capacity factor and credit distribution
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

if strcmp(stationName,'Israel wind'); figure(10008); end
hold on;
[h_ax,h_l] = plot2ylinscales(1:12,CFtempMonthly,1/CF);
 set(h_l,'Color','k'); set(h_l,'LineStyle','-');set(h_l,'lineWidth',3);
plot(1:12,CF50Monthly,'r','lineWidth',3);
plot(1:12,CF30Monthly,'b','lineWidth',3);
plot(1:12,CF1Monthly,'c','lineWidth',3);
ylabel(h_ax(1),'CC')
ylabel(h_ax(2),'CC/CF')
title({[stationName ', interannual monthly capacity credit'],['Interannual capacity factor = ' num2str(100*CF,2) '%']});
a = legend('CF','CC50','CC30','CC1');
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
end

