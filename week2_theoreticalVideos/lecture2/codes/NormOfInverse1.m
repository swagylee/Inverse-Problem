%This script will define teh example and run compute the norm of the
%convolution matrix.

clear
close all
clc

% 1 Define the function
n = 100;
f = zeros(n,1);
f(25:75) = 1;

figure(1)
plot(f,'LineWidth',3)



%2 Define the point spread function
p = [1, 2, 3, 2, 1]';
p = p/sum(p); % normalization


%3 Use the conv command to find the convolution
c = conv(f,p,'same'); % the option 'same' is not the default one

%4 Define the matrix A and the data for the inverse problem
A = convmtx(p,n);
nu = (length(p)-1)/2;
A = A(1+nu:n+nu,:); % extract only the (n) rows from 1+nu to n + nu
m = A*f;

% Naively invert this.
f_inv = inv(A)*m;

%5 Plot original sigmal and the deconvolved signal

figure(1)
plot(f,'r','LineWidth',3)
hold on
plot(f_inv,'b--','Linewidth',3)


% Add noise to Af

delta = 1E-4;    % Very small noise.
epsilon = delta*randn(n,1); % (pseudo)random gaussian numbers with mean 0 and std deviation delta 
m_d = m + epsilon;
f_inv1 = inv(A)*m_d;

figure(2)
plot(f,'r','LineWidth',3)
hold on
plot(f_inv1,'b--', 'LineWidth',3)




delta = 1E-2;
epsilon = delta*randn(n,1);
m_d1 = m + epsilon;

f_inv2 = inv(A)*m_d1;
figure(3)
plot(f,'k-','LineWidth',3)
hold on
plot(f_inv2,'b--','LineWidth',3)

delta = 1E-1;
epsilon = delta*randn(n,1);
m_d2 = m + epsilon;
f_inv3 = inv(A)*m_d2;

figure(4)
plot(f,'k-','LineWidth',3)
hold on
plot(f_inv3,'b--','LineWidth',3)

inv_A = inv(A);

norm_inv_A = norm(inv_A);

[U,S,V] = svd(A);

sing_val = diag(S);
