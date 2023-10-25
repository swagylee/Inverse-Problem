%% Exercise 1
close all
clear
clc

% 1 Create signal
n = 100;
x = (1:n)';
n0 = floor(n/15);
n1 = 2*n0; n2 = 3*n0; n3 = 4*n0; n4 = 8*n0; n5 = 10*n0; n6 = 14*n0;
f = (x>n1).*(x<=n2) + 3*(x-n3)/(n4-n3).*(x>n3).*(x<=n4)...
    - (1+cos(x*2*pi/(n6-n5))).*(x>n5).*(x<=n6);

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

% 2 Generate data
m = A*f;
delta = 1e-2;
noise = randn(n,1);
noise = noise/norm(noise);
md = m + delta*noise;

% 3 Tikhonov
alpha = 0.1;

% Strategy I: SVD
tic
[U,D,V]=svd(A);
sval = diag(D);
Ta1 = V*diag(sval./(sval.^2+alpha))*U';
fd1 = Ta1*md;
t1 = toc;

% Strategy II: normal equations
tic
fd2 = (A'*A + alpha*eye(n))\(A'*md);
t2 = toc;

% Strategy III: stacked form
tic
fd3 = [A; sqrt(alpha)*eye(n)]\([md; zeros(n,1)]);
t3 = toc;

% 4 Comparison
disp(norm(fd1-fd2))
disp(norm(fd1-fd3))
fprintf('Singular values : %.5f s \n',t1);
fprintf('Normal equations: %.5f s \n',t2);
fprintf('Stacked formula : %.5f s \n',t3);

%% Exercise 2
close all
clear
clc

% 1 Create signal
n = 200;
x = (1:n)';
n0 = floor(n/15);
n1 = 2*n0; n4 = 8*n0; n5 = 10*n0; n6 = 14*n0;
f = 5 + ((x>n1).*(x<=n4) - (1+cos(x*2*pi/(n6-n5))).*(x>n5).*(x<=n6));

% 2 Point spread function
PSF = [1 1 1];
PSF = conv(PSF,PSF);
PSF = conv(PSF,PSF);
PSF = conv(PSF,PSF);
PSF = PSF/(sum(PSF));

% Convolution matrix
A = convmtx(PSF,n);
nu = (length(PSF)-1)/2;
A = A(:,nu+1:nu+n);

% 3 Generate data
m = A*f;
delta = 1e-2;
noise = randn(n,1);
noise = noise/norm(noise);
md = m + delta*noise;

% 4 Tikhonov
alpha = 0.1;
fd = (A'*A + alpha*eye(n))\(A'*md);
figure(1)
plot(x,f,'k')
hold on
plot(x,fd,'b','LineWidth',2.5)
axis([0 n 0 7])
title('Tikhonov')

% 5 Generalized Tikhonov
f0 = 5*ones(size(f));
fdG = (A'*A + alpha*eye(n))\(A'*md + alpha*f0);
figure(2)
plot(x,f,'k')
hold on
plot(x,fdG,'g','LineWidth',2.5)
axis([0 n 0 7])
title('Generalized Tikhonov')

%% Exerise 3
close all
clear
clc

% 1 Create signal
n = 100;
x = (1:n)';
n0 = floor(n/15);
n1 = 2*n0; n4 = 8*n0; n5 = 10*n0; n6 = 14*n0;
f = (x>n1).*(x<=n4) - (1+cos(x*2*pi/(n6-n5))).*(x>n5).*(x<=n6);

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
delta = 1e-2;
noise = randn(n,1);
noise = noise/norm(noise);
md = m + delta*noise;

% 2 Tikhonov
alpha = 0.1;
fd = (A'*A + alpha*eye(n))\(A'*md);
figure(1)
plot(x,f,'k')
hold on
plot(x,fd,'b','LineWidth',2.5)
title('Tikhonov')

% 3 Generalized Tikhonov
L = eye(n) - diag(ones(1,n-1),-1); % the matrix representing the derivative
fdG = (A'*A + alpha*(L'*L))\(A'*md);
figure(2)
plot(x,f,'k')
hold on
plot(x,fdG,'g','LineWidth',2.5)
title('Generalized Tikhonov')

%% Exercise 4
close all
clear
clc

% 1 Create the signal
n = 100;
x = (1:n)';
n0 = floor(n/15);
n1 = 2*n0; n2 = 3*n0; n3 = 4*n0; n4 = 8*n0; n5 = 10*n0; n6 = 14*n0;
f = (x>n1).*(x<=n2) + 3*(x-n3)/(n4-n3).*(x>n3).*(x<=n4) - (1+cos(x*2*pi/(n6-n5))).*(x>n5).*(x<=n6);

% Point spread function from a continuous function
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
delta = 1e-3;
noise = randn(n,1);
noise = noise/norm(noise);
md = m + delta*noise;

% 2 create the sample test for alpha
ALPHA = logspace(-9,-3,50);

% 3 If I knew the solution...
error = zeros(size(ALPHA));
for i = 1:length(ALPHA)
    alpha = ALPHA(i);
    fd = (A'*A + alpha*eye(n))\(A'*md);
    error(i) = norm(fd-f)/norm(f);
end
figure(1)
loglog(ALPHA,error,'bo-');
[~,i]=min(error);
hold on
loglog(ALPHA(i),error(i),'r.','MarkerSize',25);
xlabel('\alpha')
ylabel('|| T_\alpha m_\delta - f ||')
title(['Real best one: ',num2str(ALPHA(i))])

% 4 Morozov
res = zeros(size(ALPHA));
for i = 1:length(ALPHA)
    alpha = ALPHA(i);
    fd = (A'*A + alpha*eye(n))\(A'*md);
    res(i) = norm(A*fd-md);
end
figure(2)
loglog(ALPHA,res,'bo-');
[~,i]=min(abs(res-delta));
title(['Morozov: ',num2str(ALPHA(i))])
hold on
loglog(ALPHA,delta*ones(size(ALPHA)),'--k')
loglog(ALPHA(i),res(i),'r.','MarkerSize',25);
xlabel('\alpha')
ylabel('r(\alpha)')

% 5 L curve
X = zeros(size(ALPHA));
Y = X;
for i = 1:length(ALPHA)
    alpha = ALPHA(i);
    fd = (A'*A + alpha*eye(n))\(A'*md);
    X(i) = log(norm(A*fd-md));
    Y(i) = log(norm(fd));
end
figure(3)
plot(X,Y,'bo-');

Xn=(X-min(X))/(max(X)-min(X)); %rescale X so it ranges from 0 to 1
Yn=(Y-min(Y))/(max(Y)-min(Y)); %rescale X so it ranges from 0 to 1
[~,i]=min(Xn.^2+Yn.^2); % find the closest point to the origin

title(['L curve: ',num2str(ALPHA(i))])
hold on
plot(X(i),Y(i),'r.','MarkerSize',25);
xlabel('log(r(\alpha))')
ylabel('log(|| T_\alpha m_\delta ||)')

