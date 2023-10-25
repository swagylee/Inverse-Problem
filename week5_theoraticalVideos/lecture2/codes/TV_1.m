% Total Variation example 1

% Dimension of the unknown vector 
n = 100;

% Construct the signal
f = zeros(n,1);
f(25:75) = 1;


% Point Spread function

PSF = [1 1 1];
for ii = 1:3
    PSF = conv2(PSF,PSF);
end

% Normalize PSF
PSF = PSF/sum(PSF);

figure(111)
clf
plot(f,'k','LineWidth',1.5)

figure(222)
clf
plot(PSF,'r','LineWidth',1.5)



% Construct the convolution matrix
A = convmtx(PSF,n);

% Resize the matrix
A = A(:,((length(PSF)-1)/2+1):(end-(length(PSF)-1)/2));

% Construct convolved signal that includes noise
m = A*f;

mn1 = m + 0.01*randn(size(m));
mn2 = m + 0.005*randn(size(m));
mn3 = m + 0.0001*randn(size(m));

figure(3)
clf
subplot(141)
plot(f,'k','LineWidth',1.5)
hold on
plot(m,'r','LineWidth',1.5)
subplot(142)
plot(f,'k','LineWidth',1.5)
hold on
plot(mn1,'r','LineWidth',1.5)
subplot(143)
plot(f,'k','LineWidth',1.5)
hold on
plot(mn2,'r','LineWidth',1.5)
subplot(144)
plot(f,'k','LineWidth',1.5)
hold on
plot(mn3,'r','LineWidth',1.5)










% Reconstrcut with L1 regularization
alpha = 1.7;

MAXITER = 800;
% L1_rec1 = sparsity_recon(A,eye(n),mn1,alpha,MAXITER);
% L1_rec2 = sparsity_recon(A,eye(n),mn2,alpha,MAXITER);
% L1_rec3 = sparsity_recon(A,eye(n),mn3,alpha,MAXITER);



L = -eye(n)+diag(ones(1,n-1),1);
L(end,1) = 1;

L1_rec1_d = sparsity_recon(A,L,mn1,alpha,MAXITER);
L1_rec2_d = sparsity_recon(A,L,mn2,alpha,MAXITER);
L1_rec3_d = sparsity_recon(A,L,mn3,alpha,MAXITER);



% Plot the graph

% figure(1)
% clf
% subplot(131)
% plot(f,'k','linewidth',1.2)
% hold on
% plot(mn1,'r','linewidth',1.2)
% plot(L1_rec1,'b','linewidth',1.2)
% axis square
% axis([1 n -.05 1.1]) 
% subplot(132)
% plot(f,'k','linewidth',1.2)
% hold on
% plot(mn2,'r','linewidth',1.2)
% plot(L1_rec2,'b','linewidth',1.2)
% axis square
% set(gca,'yticklabel',{})
% axis([1 n -.05 1.1]) 
% subplot(133)
% plot(f,'k','linewidth',1.2)
% hold on
% plot(mn3,'r','linewidth',1.2)
% plot(L1_rec3,'b','linewidth',1.2)
% axis square
% set(gca,'yticklabel',{})
% axis([1 n -.05 1.1]) 

figure(2)
clf
subplot(131)
plot(f,'k','linewidth',2.2)
hold on
plot(mn1,'r','linewidth',1.2)
plot(L1_rec1_d,'b','linewidth',1.2)
axis square
axis([1 n -.05 1.1]) 
subplot(132)
plot(f,'k','linewidth',2.2)
hold on
plot(mn2,'r','linewidth',1.2)
plot(L1_rec2_d,'b','linewidth',1.2)
axis square
set(gca,'yticklabel',{})
axis([1 n -.05 1.1]) 
subplot(133)
plot(f,'k','linewidth',2.2)
hold on
plot(mn3,'r','linewidth',1.2)
plot(L1_rec3_d,'b','linewidth',1.2)
axis square
set(gca,'yticklabel',{})
axis([1 n -.05 1.1]) 

