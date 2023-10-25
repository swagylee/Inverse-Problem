close all
clear
clc

n = 128;

x = linspace(0,n,n).';
f = zeros(n,1);
n0 = round(n/15);
n1 = 2*n0; n2 = 3*n0;
n3 = 4*n0; n4 = 8*n0;
n5 = 10*n0; n6 = 14*n0;
f(n1:n2) = 1;
f(n3:n4) = 3*((n3:n4) - n3)/(n4-n3);
f(n5:n6) = -1-cos((n5:n6)*2*pi/(n6-n5));


% Point spread function
PSF = [1 1 1];
PSF = conv(PSF,PSF);
PSF = conv(PSF,PSF);
PSF = conv(PSF,PSF);
PSF = PSF/(sum(PSF));

% Create convolution matrix
A = convmtx(PSF,n);
nu = (length(PSF)-1)/2;
A = A(:,nu+1:nu+n);

% Generate data
m = A*f;
delta = 0.2;
noise = randn(n,1);
md = m + delta*noise;


%% 1 ISTA for canonical basis sparsity - analysis

alpha = 0.1;
D = eye(n);
soft_thr = @(a,x) (x-a).*(x>=a) + (x+a).*(x<=-a);
r1 = ones(n,1);
for i=1:1000
    r1 = D*soft_thr(alpha,D'*(r1 + A'*(md - A*r1)));
end
figure(2)
clf
plot(1:n,r1,'b','LineWidth',2)
hold on
plot(1:n,f,'k')
plot(1:n,md,'r')
title(['L^1 norm: ', num2str(norm(f-r1)/norm(f))])

%% 2 ISTA for Haar wavelets sparsity - analysis

alpha = 0.1;
D = haarmtx(n);
soft_thr = @(a,x) (x-a).*(x>=a) + (x+a).*(x<=-a);
r2 = ones(n,1);
for i=1:1000
    r2 = D*soft_thr(alpha,D'*(r2 + A'*(md - A*r2)));
end
figure(3)
clf
plot(1:n,r2,'b','LineWidth',2)
hold on
plot(1:n,f,'k')
plot(1:n,md,'r')
title(['Haar: ', num2str(norm(f-r2)/norm(f))])

%% ISTA for a tight frame - analysis

alpha = 0.1;
D = [ haarmtx(n) eye(n)]/sqrt(2);
soft_thr = @(a,x) (x-a).*(x>=a) + (x+a).*(x<=-a);
r3 = zeros(n,1);
for i=1:1000
    r3 = D*soft_thr(alpha,D'*(r3 + A'*(md - A*r3)));
end
figure(4)
clf
plot(1:n,r3,'b','LineWidth',2)
hold on
plot(1:n,f,'k')
plot(1:n,md,'r')
title(['Tight frame - Analysis: ', num2str(norm(f-r3)/norm(f))])

%% 4 ISTA for a tight frame - synthesis

alpha = 0.05;
D = [ haarmtx(n) eye(n)]/sqrt(2);
soft_thr = @(a,x) (x-a).*(x>=a) + (x+a).*(x<=-a);
w4 = zeros(size(D,2),1);
for i=1:1000
    w4 = soft_thr(alpha,w4 + D.'*A'*(md - A*D*w4));
end
r4 = D*w4;
figure(5)
clf
plot(1:n,r4,'b','LineWidth',2)
hold on
plot(1:n,f,'k')
plot(1:n,md,'r')
title(['Tight frame - Synthesis: ', num2str(norm(f-r4)/norm(f))])

