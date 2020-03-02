% Neimark sacker bifurcation point given an alpha that minimizes a11=a33

clear; close all;

%parameters
r1 = 0.18; r2 = 0.1045;
k1 = 5e6;  k2 = 3e6; 
a2 = 3.422e-9;  
d1 = 0.0412;  d2 = 0.0412;

% Values of a1 and b from previous papers
beta_sakar = 4.32e-8; beta_kartal = 2.94e-7;a1_kartal = 2.2772e-7;

%Setup beta as free parameter
syms b
assume (b, {'positive','real'});


%variable replacements and beta limit
K1 = 1/k1; K2 = 1/k2;
b_lower = (d1*K1*K2*r2+K2*r2*a2)/(K1*(r2-d2));
%alpha simplification

a1_bar = b*r1*(K1*(-b*r1+ d1*K2*r2)+K2*r2*a2)/...
    (d2*(b*K1*r1-K2*r2*a2)+r2*((K1*(-b+d1*K2)*r1)+K2*r2*a2));

a1_prime = diff(a1_bar,b);
B_vals = double(solve(a1_prime==0,b));

b_cr = B_vals(B_vals>b_lower);

a1 = double(subs(a1_bar,b,b_cr)); 

%stable points
A = [r1*K1 a1 0; -a2 0 b; 0 b r2*K2];
c = [r1 d1 r2-d2]';
[x_interior,R] = linsolve(A,c);

vals = num2cell(x_interior);
[M,N,Z] = deal(vals{:});


E1 = exp(-r1*K1*M);
E2 = exp(-r2*K2*Z);

p2 = -1-E1-E2;

p1 = E1 + E2 + E1*E2 + b^2*N*(1- E2)/(r2*K2) - a1*a2*N*(1-E1)/(r1*K1);

p0 = -b^2*N*E1*(1-E2)/(K2*r2) - E1*E2 + a1*a2*N*E2*(1- E1)/(K1*r1);

D (b)= 1- p1 -p0^2 + p0*p2;

Beta_bifurcation = double(vpasolve(D==0,b,beta_kartal));

fprintf('The beta value from Sakar & Banerjee: %e \n',beta_sakar);
fprintf('The beta value from Kartal: %e \n',beta_kartal);
fprintf('Our calculated bifurcation beta: %e \n',Beta_bifurcation);

