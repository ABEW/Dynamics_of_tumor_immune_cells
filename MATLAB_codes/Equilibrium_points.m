% Represent the model for the tumor, hunting and resting cell interaction

syms r1 r2 k1 k2 a1 a2 b d1 d2

assume ([r1 r2 k1 k2 a1 a2 b d1 d2], {'real','positive'})
syms M N Z
assume ([M,N,Z], {'real','positive'})

g = k2*(r2-d2)/r2;

t = num2cell(zeros(9,1));
[m0,m2,m4,n0,n1,n2,n3,z0,z1] = deal(t{:});

t1 = num2cell([k1,k1,g,g,(b*r2*g-r2*d1)/(k2*b^2), d1/b]);

[m1,m3,z2,z3,n4,z4] = deal(t1{:});

Eq_pts = [m0 n0 z0; m1 n1 z1; m2 n2 z2; m3 n3 z3; m4 n4 z4];

A = [r1/k1 a1 0; -a2 0 b; 0 b r2/k2];

c = [r1 d1 r2-d2]';

[x_interior,R] = linsolve(A,c);

Var_mat (M,N,Z) = [ r1-2*r1*M/k1-a1*N -a1*M 0; -a2*N b*Z-a2*M-d1 b*N; 0 -b*Z ...
    r2-2*r2*Z/k2-b*N-d2];

Eig_vals =[];

for i=1:size(Eq_pts,1)
    vals = num2cell(Eq_pts(i,:));
    [x,y,z] = deal(vals{:});
    Eig_vals = [Eig_vals eig(Var_mat(x,y,z))];
end

vals = num2cell(x_interior);
[x,y,z] = deal(vals{:});
Eig_vals = [Eig_vals eig(Var_mat(x,y,z))];
