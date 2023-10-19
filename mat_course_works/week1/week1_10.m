p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = p/sum(p); 

g = @(x) (2<x & x<=5).*(1)+(5<x & x<=8).*(0.5)+...
    (8<x & x<=11).*(0.75)+(11<x & x<=14).*(0.25);

x1 = linspace(0,16,1000);
x2 = linspace(0,16,32);
x3 = linspace(0,16,16);

f1 = g(x1)';
f2 = g(x2)';
f3 = g(x3)';

M = (length(p)-1)/2;

A1 = convmtx(p,length(f1));
A1 = A1((M+1):(end-M),:);

A2 = convmtx(p,length(f2));
A2 = A2((M+1):(end-M),:);

A3 = convmtx(p,length(f3));
A3 = A3((M+1):(end-M),:);

rng(0,'twister');
noise1 = 1/norm(f1)*randn(size(A1*f1));
m0_1 = A1*f1;
m1_1 = A1*f1+0.01*noise1;
m2_1 = A1*f1+0.1*noise1;

noise2 = 1/norm(f2)*randn(size(A2*f2));
m0_2 = A2*f2;
m1_2 = A2*f2+0.01*noise2;
m2_2 = A2*f2+0.1*noise2;

noise3 = 1/norm(f3)*randn(size(A3*f3));
m0_3 = A3*f3;
m1_3 = A3*f3+0.01*noise3;
m2_3 = A3*f3+0.1*noise3;
% Perform the deconvolution
f1_0 = inv(A1)*m0_1;
f1_1 = inv(A1)*m1_1;
f1_2 = inv(A1)*m2_1;

% Perform the deconvolution
f2_0 = inv(A2)*m0_2;
f2_1 = inv(A2)*m1_2;
f2_2 = inv(A2)*m2_2;

% Perform the deconvolution
f3_0 = inv(A3)*m0_3;
f3_1 = inv(A3)*m1_3;
f3_2 = inv(A3)*m2_3;
% Calculate the relative errors
relative_error1_0 = norm(f1 - f1_0) / norm(f1) * 100
relative_error1_1 = norm(f1 - f1_1) / norm(f1) * 100
relative_error1_2 = norm(f1 - f1_2) / norm(f1) * 100

relative_error2_0 = norm(f2 - f2_0) / norm(f2) * 100
relative_error2_1 = norm(f2 - f2_1) / norm(f2) * 100
relative_error2_2 = norm(f2 - f2_2) / norm(f2) * 100

relative_error3_0 = norm(f3 - f3_0) / norm(f3) * 100
relative_error3_1 = norm(f3 - f3_1) / norm(f3) * 100
relative_error3_2 = norm(f3 - f3_2) / norm(f3) * 100
