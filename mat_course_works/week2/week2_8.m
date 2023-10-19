p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = p / sum(p);

M = (length(p)-1)/2;

n = 1024;

A = convmtx(p, n);

A = A((M+1):(end-M),:);

x1 = linspace(0, 16, n);
g = @(x) (2<x & x<=5).*(1) + (5<x & x<=8).*(0.5) + (8<x & x<=11).*(0.75) + (11<x & x<=14).*(0.25);
f = g(x1)';

rng(0, 'twister');
noise = 1/norm(f) * randn(size(A * f));
m0 = A * f;
m1 = A * f + 0.00001 * noise;
m2 = A * f + 0.01 * noise;

B = pinv(A);
finvA_0 = inv(A) * m0;
finvA_1 = inv(A) * m1;
finvA_2 = inv(A) * m2;

fB_0 = B * m0;
fB_1 = B * m1;
fB_2 = B * m2;

norm_diff_0 = norm(finvA_0 - fB_0)
norm_diff_1 = norm(finvA_1 - fB_1)
norm_diff_2 = norm(finvA_2 - fB_2)


