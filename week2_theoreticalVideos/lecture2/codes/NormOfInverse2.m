clear
close all
clc

%1 The function f
n = 100;
f = zeros(n,1);
f(25:75) = 1;

% Point spread function
p = [1, 1, 1, 1, 1]';
p = conv(p,p); % full autoconvolution
p = p/sum(p);

%2 Find the convolution matrix A
A = convmtx(p,n);
nu = (length(p)-1)/2;
A = A(nu+1:nu+n,:);

%3 Data for the inverse problem
m = A*f;
f0 = inv(A)*m;


figure(1)
plot(m,'k-','LineWidth',3)
hold on
plot(f0,'g--','LineWidth',3)


%4 Add noise to the data.
delta = 1E-4;
epsilon = delta*randn(n,1); % (pseudo)random gaussian numbers with mean 0 and std deviation delta 
m_d = m + epsilon;
f_d = inv(A)*m_d;
figure(2)
plot(f,'k-','LineWidth',3)
hold on
plot(f_d,'g--','LineWidth',3)

%4 More noise
delta = 1E-2;
epsilon = delta*randn(n,1);
m_d1 = m + epsilon;
f_d1 = inv(A)*m_d1;
figure(3)
plot(f,'k-','LineWidth',3)
hold on
plot(f_d1,'g--','LineWidth',3)
error_d1 = norm(f-f_d1);

%5 Gaussian distribution of the point spread function
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
f_d2 = inv(A2) * m_d2;
figure(4)
plot(f,'k-','LineWidth',3)
hold on
plot(f_d2,'g--','LineWidth',3)



inv_A = inv(A);

norm_inv_A = norm(inv_A);

[U,S,V] = svd(A);

sing_val = diag(S);

cond_A = cond(A);
