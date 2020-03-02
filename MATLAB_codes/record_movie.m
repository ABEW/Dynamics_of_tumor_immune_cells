clear; close all;
init = [2e5,2e5,1e5];
b = 3e-8;
b_scaled = b*1e8;

f = figure;
ax = axes('Parent',f,'position',[0.13 0.2  0.77 0.7]);
start = plot_helper(init,b_scaled); 
P = plot3(start(:,1),start(:,2),start(:,3));
grid('on');
xlabel('Hunting cell population (N)','FontSize',15);
ylabel('Tumor cell population (M)','FontSize',15);
zlabel('Resting cell population (Z)','FontSize',15);
title(['\beta = ',num2str(b)],'FontSize',15);
set(get(gca,'xlabel'),'rotation',10);
set(get(gca,'ylabel'),'rotation',-15);
    
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
xlim([0,inf]);
ylim([0,inf]);
zlim([0,inf]);

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

loops = 200;
F(loops) = struct('cdata',[],'colormap',[]);
b_val = linspace(4e-8,6e-7,loops);

for j = 1:loops
    update_figure(ax,P,[2e5,2e5,1e5],b_val(j));
    F(j) = getframe(gcf);
end

v = VideoWriter('bifurcation','Uncompressed AVI');
open(v);
writeVideo(v,F);
close(v);

% fig = figure;
% movie(fig,F)

