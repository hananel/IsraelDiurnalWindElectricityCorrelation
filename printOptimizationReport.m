
figList = get(0,'Children');
for i=figList'
    figure(i);
    saveas(gcf,[reportDir 'fig',num2str(i),'.fig'])
    print('-dpsc','-append',[reportDir reportName])
end