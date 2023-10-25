% We experiment with the naive inversion, SVD and inversion with TSVD.




% Choose the dimension of the unknown vector 
N = 62;

% Construct the evealuation points for the signal.
x = linspace(0,2*pi,N);

% Construct the signal
f = (2*pi - x) .* x.^3 .*(1-cos(x));

% Plot this signal to get a look

figure(1)
plot(x,f,'k','LineWidth',2)


% Construct the point spread function

p = [1 1 1];

for ii = 1:3
    p = conv2(p,p);  % Autoconvolution
end

p = p/sum(p);
% Plot the Point Spread Function

figure(2)
plot(p,'r')
hold on 
plot(p,'r*','markersize',10)


% Construct the convolution matrix

A = convmtx(p,N);

% Remove the extra columns of the convolved matrix
nu = (length(p)-1)/2;

A = A(:,nu+1:end-nu);

% Compute the SVD of matrix A
[U,D,V] = svd(A);
svals = diag(D);

% Lets look at the matrix A

figure(3)
subplot(121)
spy(A)
title(['Convolution Matrix of size',num2str(size(A,1)),'x',num2str(size(A,2))])
subplot(122)
semilogy(svals)
title(['Semilog plot of singular values'])


% Compute the convolved signal that includes noise.

% First signal without noise
m = A*f(:);

% Signal with noise
mn = A*f(:) + randn(size(m));

% Plot the original signal and the convolved signal

figure(4)
clf 
plot(x,f,'k','LineWidth',1.5)
hold on
%plot(x,m,'b','LineWidth',1.8)
plot(x,mn,'r','LineWidth',1.2)



% Naive inversion

f_0 = inv(A)*mn;

% Using pseudoinverse.b

f_1 = pinv(A)*mn;

% Plot the results

figure(5)
clf
plot(x,f,'k','LineWidth',1.5)
hold on
plot(x,mn,'r','LineWidth',1.2)
plot(x,f_0,'b','LineWidth',2)
plot(x,f_1,'g','LineWidth',2.5)



% Lets truncate the singular values and do reconstruction

% Choose the value of r_alpha
D_ralpha = zeros(size(D.'));
r_alpha = 10;
dralpha = svals(1:r_alpha);
D_ralpha(1:r_alpha,1:r_alpha) = diag(1./dralpha);





% Reconstruct with TSVD

f_TSVD = V*D_ralpha*U.'*mn;

% Use pinv(A, 'tol')
f_TSVD1 = pinv(A,1E-2)*mn;
figure(6)
clf
subplot(121)
plot(x,mn,'k','LineWidth',1.5)
hold on
plot(x,f_TSVD,'b','LineWidth',2)
%plot(x,f_TSVD1,'g','LineWidth',1.8)
subplot(122)
plot(V(:,10),'k')
