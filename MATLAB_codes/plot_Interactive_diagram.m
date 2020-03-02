function [P] = plot_Interactive_diagram(init, b)
% test ploting for time delayed de

r1 = 0.18; r2 = 0.1045;
k1 = 5e6;  k2 = 3e6; 
a2 = 3.422e-9; %b = 4.32e-8; a1 = 2.2772e-7;
d1 = 0.0412;  d2 = 0.0412;

% Swap K values and set beta limits
K1 = 1/k1; K2 = 1/k2;

b_Upper = K2*r2*(K1+K1*d1+a2)/(K1*(r2-d2));
b_Lower = (K1*K2*r2*d1+K2*r2*a2)/(K1*(r2-d2));
b_scaled = 1e8*b;


slider_max = round(1e8*b_Upper,2);
slider_min = round(1e8*b_Lower*0.8,2);


% Setup figure and slider control
f = figure;
ax = axes('Parent',f,'position',[0.13 0.3  0.77 0.64]);
start = plot_helper(init,b_scaled); 
P = plot3(start(:,1),start(:,2),start(:,3));
grid('on');
xlabel('Hunting cell population (N)');
ylabel('Tumor cell population (M)');
zlabel('Resting cell population (Z)');
title(['\beta = ',num2str(b)]);

h = uicontrol('Parent',f,'Style','slider','Position',[81,40,419,23],...
              'value',b_scaled, 'min',slider_min, 'max',slider_max);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,40,25,25],...
                'String',slider_min,'BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,40,35,25],...
                'String',slider_max,'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','Conversion rate (beta)','BackgroundColor',bgcolor);

%txh = uicontrol('Parent',f,'Style','text','Position',[275,10,50,15],...
   % 'BackgroundColor','w');
%set(txh,'String',num2str(round(b_scaled),2));

h.Callback = @(es,ed) updateSystem(es.Value,ax,P,plot_helper(init,es.Value));

end