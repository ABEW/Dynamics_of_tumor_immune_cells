function[] = update_figure(ax,P,init, b)
% Takes initial values for M,N and Z and generates the plot 

r1 = 0.18; r2 = 0.1045;
k1 = 5e6;  k2 = 3e6; 
a2 = 3.422e-9; 
d1 = 0.0412;  d2 = 0.0412;

%variable replacements
K1 = 1/k1; K2 = 1/k2;

% WHEN a11 = a33
%       a1 = b*r1*(K1*(-b*r1+ d1*K2*r2)+K2*r2*a2)/...
%       (d2*(b*K1*r1-K2*r2*a2)+r2*((K1*(-b+d1*K2)*r1)+K2*r2*a2));

a1 = 2.2683e-7;

M = init(1); N = init(2); Z = init(3);

for n=1:2000
    c1 = (r1-a1*N(n));
    c2 = r1*K1;
    c3 = r2-b*N(n)-d2;
    c4 = r2*K2;
    
    m_n1 = M(n)*c1/((c1-c2*M(n))*exp(-c1)+c2*M(n));
    n_n1 = N(n) * exp(b*Z(n)-d1-a2*M(n));
    z_n1 = Z(n)*c3/((c3-c4*Z(n))*exp(-c3)+c4*Z(n));
    
    M = [M, m_n1];
    N = [N, n_n1];
    Z = [Z, z_n1];
end
set(P,'xdata',N);
set(P,'ydata',M);
set(P,'zdata',Z);
title(ax,['\beta = ',num2str(b,3)]);


end