function [] = updateSystem(b_scaled,ax,P,V)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
set(P,'xdata',V(:,1));
set(P,'ydata',V(:,2));
set(P,'zdata',V(:,3));
title(ax,['\beta = ',num2str(round(b_scaled,2)*1e-8,3)]);
%set(txh,'String',num2str(beta))
drawnow;
end

