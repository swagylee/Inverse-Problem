% Here we consider the discrete convolution matrix given by f[25:75] = 1
% and zero otherwise, where n = [1:100]

% We consider zeropadding point spread function given by p =
% [1,1,1,1,1]' and normalize it. 



%1 The function f
n = 100;
f = zeros(n,1);
f(25:75) = 1;

% Point spread function
p = [1, 1, 1, 1, 1]';
p = conv(p,p); % full autoconvolution
p = p/sum(p);  % Normalize the point spread function.

%2 Find the convolution matrix A
A = convmtx(p,n);
nu = (length(p)-1)/2;
A = A(nu+1:nu+n,:);

% Let's create blurred image/signal given m with inverse crime

m = A*f;
m_n = A*f +0.01* randn(n,1);
% Plot the original signal f and blurred signal m

% figure(1)
% plot(f,'k','LineWidth',1.5)
% hold on
% plot(m,'b','LineWidth',1.5)

% To find the original signal given m,

f_0 = inv(A)*m;

f_1 = pinv(A) *m;

% figure(2)
% plot(f,'k','LineWidth',1.5)
% hold on
% plot(f_0,'b','LineWidth',1.5)
% plot(f_1,'r','LineWidth',1.5)

% Add noise to the signal

delta = 1E-4;
epsilon = delta * randn(n,1);

m_d = A*f + epsilon;

f_noise_0 = inv(A) * m_d;
f_noise_1 = pinv(A) *m_d;

 figure(3)
 subplot(121)
 plot(f,'k','LineWidth',1.5)
 hold on
 plot(f_noise_0,'b','LineWidth',1.5)
 subplot(122)
 plot(f,'k','LineWidth',1.5)
 hold on
 plot(f_noise_1,'r','LineWidth',1.5)
 
 % Analyse the matrix
 
 [U D V] = svd(A);
 d = diag(D);
 
 % Plotting of the analysis of the matrix
 
 figure(4)
 clf 
 subplot(121)
 spy(A)
 subplot(122)
 semilogy(d)
 
 f_noise_3 = pinv(A,1E-1)*m;
 
 
 figure(5)
 clf
 plot(f,'k','LineWidth',2)
 hold on
 plot(f_noise_3,'r','LineWidth',1.5)
 
 figure(6)
 clf 
 plot(V(:,100))
 
%  ralpha = 30;
%  d_ralpha = diag(1:ralpha);
%  D_plus = zeros(size(D.'));   % lets define a zero matrix of size D
%  D_plus(1:ralpha,1:ralpha) = diag(1./d_ralpha);
%  
%  
%  
%  % Let's anaylse the matrix A
%   [U,D,V] = svd(A);
%   d = diag(D);
%   
%   
%  figure(4)
%  subplot(121)
%  spy(A)
%  subplot(122)
%  semilogy(d)
%  
%  f_noise_2 = pinv(A,1E-3)*m_d;
%  
%  figure(5)
%  clf
%  subplot(121)
%  plot(f,'k','LineWidth',1.5)
%  hold on
%  plot(f_noise_1,'b','LineWidth',1.5)
%  subplot(122)
%  plot(f,'k','LineWidth',1.5)
%  hold on
%  plot(f_noise_2,'r','LineWidth',1.5)
%  
%  
 D = zeros(size(D.'));   % Define the matrix D_plus to have the same size as that of D
 ralpha = 50;         % Number of singular values to be taken into account
 dalpha = d(1:ralpha); % Only the first ralpha diagonal values are considered.
 D(1:ralpha,1:ralpha) = diag(1./dalpha);
 f_noise_4 = V*D*U.' *m_d;
 
 
 % Plot of singular vectors
 
 figure(60)
 clf 
 subplot(121)
 plot(V(:,50))
 %subplot(122)
 subplot(122)
 plot(f,'k','LineWidth',1.5)
 hold on
 plot(f_noise_4,'r','LineWidth',1.5)
  
