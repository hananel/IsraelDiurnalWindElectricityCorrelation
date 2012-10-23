function plotStationLocations(stationVec,col,weights)

meta = loadMeta;
plotIsraelBorders
axis tight; hold on;

for i=1:length(meta)
    % checking for anemometer
    yesAnemometer = not(or(strcmp(meta(i).anemometer,'No'),meta(i).h(2)==-1));
    if yesAnemometer
        plot(meta(i).long,meta(i).lat,'.','color','g');
    end

end    
counter = 0;
for i=stationVec
    counter = counter + 1;
    if nargin > 2
        % plotting weights - in percents
        if abs(weights(counter)*100)>0
            plot(meta(i).long,meta(i).lat,'o','MarkerSize',abs(weights(counter)*100),'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
            text(meta(i).long-0.05,meta(i).lat,num2str(meta(i).num),'color','k','fontweight','bold')
        end
    else
        text(meta(i).long-0.05,meta(i).lat,num2str(meta(i).num),'color',col(counter,:),'fontweight','bold')
    end
end 
xlabel('longitude [deg]'); ylabel('latitude [deg]'); title('IMS stations used in analysis');

