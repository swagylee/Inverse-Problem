%% Example 1
clear
close all
clc

% 1
n = 100;
f = zeros(n,1);
f(25:75) = 1;

figure(1)
plot(f,'LineWidth',3)

%2
p = [1, 2, 3, 2, 1]';
p = p/sum(p); % normalization

%3
c = conv(f,p,'same'); % the option 'same' is not the default one

%4
A = convmtx(p,n);
nu = (length(p)-1)/2;
A = A(1+nu:n+nu,:); % extract only the (n) raws from 1+nu to n + nu
c1 = A*f;

%5
figure(1) 
hold on
plot(c,'r','LineWidth',3)
plot(c1,'g','LineWidth',3)
error = norm(c-c1); % quantitative comparison between c and c1

%6
% the command A = diag(v) creates a matrix having the vector v on the
% main diagonal. diag(v,1) puts v on the first upper diagonal. 
% Pay attentiont to the size of upper and lower diagonals!
A2 = p(5)*diag(ones(n-2,1),-2) + ...
     p(4)*diag(ones(n-1,1),-1) + ...
     p(3)*diag(ones(n,1))      + ...
     p(2)*diag(ones(n-1,1),1)  +...
     p(1)*diag(ones(n-2,1),2);

figure(2)
spy(abs(A-A2)) % graphical comparison: spy(M) displays the non-null elements of M
errorA = norm(A-A2); % quantitative comparison (this is a matrix norm!)
disp(errorA)

%% Example 2
clear
close all
clc

%1
n = 100;
f = zeros(n,1);
f(25:75) = 1;
p = [1, 1, 1, 1, 1]';
p = conv(p,p); % full autoconvolution
p = p/sum(p);

%2
A = convmtx(p,n);
nu = (length(p)-1)/2;
A = A(nu+1:nu+n,:);

%3
m = A*f;
f0 = A\m; % solves A*f0 = m  !! we don't need inv(A) !!
figure(1)
plot(f,'k-','LineWidth',3)
hold on
plot(f0,'g--','LineWidth',3)
error0 = norm(f-f0);

%4
delta = 1E-4;
epsilon = delta*randn(n,1); % (pseudo)random gaussian numbers with mean 0 and std deviation delta 
m_d = m + epsilon;
f_d = A\m_d;
figure(2)
plot(f,'k-','LineWidth',3)
hold on
plot(f_d,'g--','LineWidth',3)
error_d = norm(f-f_d);


%4
delta = 1E-2;
epsilon = delta*randn(n,1);
m_d1 = m + epsilon;
f_d1 = A\m_d1;
figure(3)
plot(f,'k-','LineWidth',3)
hold on
plot(f_d1,'g--','LineWidth',3)
error_d1 = norm(f-f_d1);

%5
p2 = [1, 1, 1, 1, 1]';
p2 = conv(p2,p2);
p2 = conv(p2,p2);
p2 = conv(p2,p2);
p2 = p2/sum(p2);
A2 = convmtx(p2,n);
nu = (length(p2)-1)/2;
A2 = A2(nu+1:nu+n,:);
delta = 1E-4;
epsilon = delta*randn(n,1);
m_d2 = A2*f + epsilon;
f_d2 = A2\m_d2;
figure(4)
plot(f,'k-','LineWidth',3)
hold on
plot(f_d2,'g--','LineWidth',3)
error_d2 = norm(f-f_d2);









