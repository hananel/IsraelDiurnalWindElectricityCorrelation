function plotIsraelBorders()

israelBorders0 = csvread('IsraelBorders0.csv');
israelBorders1 = csvread('IsraelBorders1.csv');
israelBorders2 = csvread('IsraelBorders2.csv');

plot(israelBorders1(1:2:end),israelBorders1(2:2:end),'color',[0.1,0.1,0.1]); 
hold on;
plot(israelBorders2(1:2:end),israelBorders2(2:2:end),'color',[0.1,0.1,0.1]);
plot(israelBorders0(:,1),israelBorders0(:,2),'color',[0.1,0.1,0.1],'linewidth',2);
axis equal