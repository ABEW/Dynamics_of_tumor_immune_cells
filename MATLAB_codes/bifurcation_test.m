% Neimark sacker bifurcation point, plot beta as function of alpha

clear; close all;

%parameters
r1 = 0.18; r2 = 0.1045;
k1 = 5e6;  k2 = 3e6; 
a2 = 3.422e-9;
d1 = 0.0412;  d2 = 0.0412;

syms b
assume (b, {'positive','real'});

%variable replacements
K1 = 1/k1; K2 = 1/k2;

sample_points = 100;
a_range = linspace(1e-8,1e-6,sample_points);
valid_index =zeros(sample_points,1);
Beta = zeros(sample_points,1);

% Set the upper and lower limits on beta
beta_lower = (d1*K1*K2*r2+K2*r2*a2)/(K1*(r2-d2));  % 3.2e-8
beta_upper = K2*r2*(K1+K1*d1+a2)/(K1*(r2-d2));   % 5.8e-7

for i = 1:length(a_range)
%stable points
    a1 = a_range(i);
    
    %Solve for the equilibtrium point
    M = (b^2*r1+a1*(b*(d2-r2)+d1*K2*r2))/(b^2*K1*r1-K2*r2*a1*a2);
    N = (r1/a1)*(1-K1*M);
    Z = (d1+a2*M)/b;

    % Solve for bifurcation point with Schur-Cohn criterior
    E1 = exp(-r1*K1*M);
    E2 = exp(-r2*K2*Z);
    p2 = -1-E1-E2;
    p1 = E1 + E2 + E1*E2 + b^2*N*(1- E2)/(r2*K2) - a1*a2*N*(1-E1)/(r1*K1);
    p0 = -b^2*N*E1*(1-E2)/(K2*r2) - E1*E2 + a1*a2*N*E2*(1- E1)/(K1*r1);

    D = simplify(1- p1 -p0^2 + p0*p2);

    b_val = double(vpasolve(D==0,b));

    % Filter out significantly erroneous results 
    if abs(b_val-Beta(i-1))>5e-6
        b_val = double(vpasolve(D==0,b,Beta(i-1)));
    end
    
    if (b_val>beta_lower && b_val<beta_upper)
        valid_index(i)=1;
    end
    
    Beta(i) = b_val;

end

a_valid = a_range(valid_index); b_valid = Beta(valid_index);
figure
hold on
plot(a_valid, b_valid,'.')
yline(3.2e-8,'r','DisplayName','Beta lower bound');
yline(5.8e-7,'g','DisplayName','Beta upper bound');
legend

[a_new,b_new] = filter_bifurcation_vals(a_valid,b_valid);


