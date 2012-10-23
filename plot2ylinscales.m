function  [h_ax,h_l] = plot2ylinscales(x,y,fac)
%PLOT2YLINSCALES - XY-plot with two linearly related Y-scales
%PLOT2YLINSCALES(X,Y,FUN_STR) - plots the data in Y vs. X and adds a
% second Y-scale which is related to the original Y data by the relationship
% specified in FUN_STR.
% FUN_STR must be a string which represents a linear function,
% e.g. 'a*x + b', where the constants 'a' and 'b' must be defined beforehand.
%
%WARNING! If FUN_STR does not evaluate to a linear function
%the second scale will be WRONG!

lin_fun = inline([num2str(fac) '*x']);
y2 = lin_fun(y);

[ax,h1,h2]=plotyy(x,y,x,y2);

delete(h2)
set(ax(1),'box','off','ytickmode','auto');
set(ax(2),'ylim',lin_fun(get(ax(1),'ylim')));
set(ax(2),'ytickmode','auto');

if nargout~=0
    h_ax = ax;h_l = h1;
end
end

