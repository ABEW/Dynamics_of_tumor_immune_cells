function [VAL] = plot_helper(init,b_scaled)
%Used to update the plot when beta value is changed

r1 = 0.18; r2 = 0.1045;
k1 = 5e6;  k2 = 3e6; 
a2 = 3.422e-9; %a1 = 2.2772e-7;
d1 = 0.0412;  d2 = 0.0412;

% Swap K values and set beta limits
K1 = 1/k1; K2 = 1/k2;
b = b_scaled*1e-8;

a1 = 2.2683e-7;

M = init(1);
N = init(2);
Z = init(3);

for n=1:50000
    c1 = (r1-a1*N(n));
    c2 = r1/k1;
    c3 = r2-b*N(n)-d2;
    c4 = r2/k2;
    
    m_n1 = M(n)*c1/((c1-c2*M(n))*exp(-c1)+c2*M(n));
    n_n1 = N(n) * exp(b*Z(n)-d1-a2*M(n));
    z_n1 = Z(n)*c3/((c3-c4*Z(n))*exp(-c3)+c4*Z(n));
    
    M = [M, m_n1];
    N = [N, n_n1];
    Z = [Z, z_n1];
end

VAL = [N',M',Z'];

end

