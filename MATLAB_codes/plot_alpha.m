%clear; close all

% syms b
% assume(b,{'positive','real'})

r1 = 0.18; r2 = 0.1045;
k1 = 5e6;  k2 = 3e6; 
a2 = 3.422e-9; 
d1 = 0.0412;  d2 = 0.0412;

% change of variables 
K1 = 1/k1; K2 = 1/k2;

bif_line = fit(a_new',b_new,'linearinterp');

b = 3e-8:1e-9:6e-7;

% When a11 = a33
a1_eq = r1*(K1*(-b.^2*r1+ b*d1*K2*r2)+b*K2*r2*a2)./...
    (d2*(b*K1*r1-K2*r2*a2)+r2*((K1*(-b+d1*K2)*r1)+K2*r2*a2));

% Limits for existence of interior equilibrium
a1_upper_limit = b.^2*r1./(b*(r2-d2)-r2*d1*K2);
a1_upper_limit_2 = b.^2*r1*K1/(a2*K2*r2); % will not be plotted
beta_lower = (d1*K1*K2*r2+K2*r2*a2)/(K1*(r2-d2));  % 3.2e-8

% Local stability equirement
beta_upper = K2*r2*(K1+K1*d1+a2)/(K1*(r2-d2));   % 5.8e-7

% Beta independence of alpha limit
a1_strict_limit = 4*d1*K2*r1*r2/ (d2-r2)^2;

% Selected alpha
a1_selected = 2.2683e-7;

x_values = 0:1e-8:16e-7;

figure
hold on
plot(a1_eq,b,'b','DisplayName','a_{11}=a_{33}')
plot(a1_upper_limit,b,'r-.','DisplayName','\alpha_1 E_5 existence limit')
plot(x_values,feval(bif_line,x_values)')
a_choose = xline(a1_selected,'k:','Chosen \alpha_1');
a_strict = xline(a1_strict_limit,'g-.',{'\alpha_1 strict','upper limit'});
b_low_line = yline(beta_lower,'k-.','\beta E_5 existence limit');
b_high_line = yline(beta_upper,'k-.','\beta E_5 local stability limit');

a_choose.LabelHorizontalAlignment = 'left';
a_strict.LabelHorizontalAlignment = 'right';
b_low_line.LabelHorizontalAlignment = 'center';
b_high_line.LabelHorizontalAlignment = 'center';

ylim([0,7e-7])
legend



