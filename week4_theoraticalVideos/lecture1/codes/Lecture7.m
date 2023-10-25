close all
clear
clc

% 1 Create the signal
n = 200;
x = linspace(0,n,n).';
f = zeros(n,1);
n0 = round(n/15);
n1 = 2*n0; n2 = 3*n0;
n3 = 4*n0; n4 = 8*n0;
n5 = 10*n0; n6 = 14*n0;
f(n1:n2) = 1;
f(n3:n4) = 3*((n3:n4) - n3)/(n4-n3);
f(n5:n6) = -1-cos((n5:n6)*2*pi/(n6-n5));

figure(1)
plot(x,f,'k','LineWidth',2)
axis([0 n -2 3])
title('Original signal')

%2 Point Spread Function and convolution matrix
PSF = [1 4 8 16 19 15 10 7 1];
PSF = PSF/sum(PSF);
figure(2)
plot(PSF,'ro-')
title('Point spread function')

A = convmtx(PSF,n);
nu = (length(PSF)-1)/2;
A = A(:,(nu+1):(nu+n));

%3 Noisy measurements
delta = 1e-5;
m = A*f;
noise = randn(size(m));
noise = noise/norm(noise);
md = m + delta*noise;

figure(3)
plot(x,f,'k')
hold on
plot(x,md,'r')
title('Noisy data')
 
%4 Pseudoinverse
psA = pinv(A,0);
fd = psA*md;

figure(4)
plot(x,fd,'k')
hold on
plot(x,f,'b')

%5 Tikhonov regularization
alpha = 1e-1;
[U,D,V] = svd(A);
svals = diag(D);
Dp_a = diag(svals./(svals.^2+alpha));
T_a = V*Dp_a*U.';
fd_a = T_a*md;

figure(5)
plot(x,fd_a,'k')
hold on
plot(x,f,'b')

%6 Comparison
error_PS = norm(fd-f)/norm(f);
error_Tik = norm(fd_a-f)/norm(f);
fprintf(' Error with pseudo-inverse: \t %.3f \n Error with Tikhonov regul: \t %.3f \n',error_PS,error_Tik)

figure(6)
subplot(1,2,1)
plot(x,fd_a,'g')
hold on
plot(x,f,'b')
axis([0 n -2 3])
title(['Tikhonov error: ',num2str(error_Tik)])
subplot(1,2,2)
plot(x,fd,'k')
hold on
plot(x,f,'b')
axis([0 n -2 3])
title(['Pseudoinverse error: ',num2str(error_PS)])

% 6 
% Modify the values of alpha and delta and run the code again
