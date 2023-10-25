%% Ex 1

close all
clear 
clc

% 1 Original signal
img = imread('mri.tif');
f = im2double(img);
figure(1)
imshow(f)
truesize([400 400])
title('Original Image');

% 2 Filters
px = 7; 
py  = 7;
p_Average = fspecial('average',[px,py]);
p_Gaussian = fspecial('gaussian',[px,py],2);
p_Motion = fspecial('motion',px,0);

% 3 Convolution
c_Average = conv2(f,p_Average,'same');
figure(2)
imshow(c_Average)
truesize([400 400])
title('Average filter');

c_Gaussian = conv2(f,p_Gaussian,'same');
figure(3)
imshow(c_Gaussian)
truesize([400 400])
title('Gaussian filter');

c_Motion = conv2(f,p_Motion,'same');
figure(4)
imshow(c_Motion)
truesize([400 400])
title('Motion filter');

% 4 Noise
mu = 0;
delta = 0.001;

md_Average = imnoise(c_Average, 'gaussian', mu, delta);
figure(2)
imshow(md_Average)
truesize([400 400])
title('Average filter + noise');

md_Gaussian = imnoise(c_Gaussian, 'gaussian', mu, delta);
figure(3)
imshow(md_Gaussian)
truesize([400 400])
title('Gaussian filter + noise');

md_Motion = imnoise(c_Motion, 'gaussian', mu, delta);
figure(4)
imshow(md_Motion)
truesize([400 400])
title('Motion filter + noise');

%% Ex 2

close all
clear 
clc

% 1 Original signal
img = imread('cameraman.tif');
F = im2double(img);
figure(1)
imshow(F)
truesize([400 400])
title('Original Image');

% 2 Filter
px = 5; 
py  = 5;
P = fspecial('average',[px,py]);

% 3 Convolution via conv2
C1 = conv2(F,P,'same');
figure(2)
imshow(C1)
truesize([400 400])
title('conv2');

% 4 Convolution matrix
[n,m] = size(F);
nu = (px-1)/2;
mu = (py-1)/2;
Nr = (n+2*nu)*(m+2*mu);

A = convmtx2(P,n,m);
 
I = [1:(n+2*nu)*mu, (Nr-mu*(n+2*nu)+1):Nr];
for i = 1:nu
    I = [I, i:(n+2*nu):Nr, (n+2*nu-i+1):(n+2*nu):Nr];
end
A(I,:) = [];

% 5 Convolution via A
c2 = A*F(:);
C2 = reshape(c2,[n,m]);
figure(3)
imshow(C2)
truesize([400 400])
title('convmtx2');

% 6 Compare
figure(4)
imshow(abs(C1-C2))
truesize([400 400])
title(['Difference: ',num2str(norm(C1-C2))]);


%% Ex 3

close all
clear 
clc

% Original signal
% img = imread('cameraman.tif');
img = imread('logo.tif');
f = im2double(img);

[n,m] = size(f);
figure(1)
imshow(f);
truesize([400 400])
title('Original Image');
set(gca,'FontSize',20)

% Filter
px = 5; 
py  = 5;
nux = (px-1)/2;
nuy = (py-1)/2;
P = fspecial('gaussian',[px,py],2);

% Convolution matrix
A = convmtx2(P,n,m);
Nc = (n+px-1)*(m+py-1);
I = [1:(n+px-1)*nuy,(Nc-nuy*(n+px-1)+1):Nc];
for i = 1:nux
    I = [I, i:(n+px-1):Nc, (n+px-i):(n+px-1):Nc];
end
A(I,:) = [];

% Noisy measurements
mu = 0;
delta = 0.001;
md = imnoise(conv2(f,P,'same'), 'gaussian');
figure(2)
imshow(md);
truesize([400 400])
title('Blurred and noisy data')
set(gca,'FontSize',20)

%Naive inversion
R1 = A\md(:);
r1 = reshape(R1,[n,m]);
figure(3)
imshow(r1);
truesize([400 400])
title('No regularization')
set(gca,'FontSize',20)

%Tikhonov
alpha = 0.15;
R2 = (A.'*A + alpha*speye([n*m,n*m]))\(A.'*md(:)); %use speye instead of eye to guarantee sparsity
r2 = reshape(R2,[n,m]);
figure(4)
imshow(r2);
truesize([400 400])
title(['Tikhonov, \alpha = ',num2str(alpha)])
set(gca,'FontSize',20)

%Generalized Tikhonov
alpha = 0.19;

Lx = spdiags([ones(n*m,1) -ones(n*m,1)], [0,-n], n*m, n*m);
Ly = spdiags([ones(n*m,1) -ones(n*m,1)], [0,1], n*m, n*m);
Ly((m:m:n*m-1)+1,(m:m:n*m)) = 0;

R3 = (A'*A + alpha*(Lx'*Lx) + alpha*(Ly'*Ly))\(A'*md(:));
r3 = reshape(R3,[n,m]);
figure(5)
imshow(r3);
truesize([400 400])
title(['Generalized Tikhonov, \alpha = ',num2str(alpha)])
set(gca,'FontSize',20)


%% WARNING: this is very slow on 256x256 pictures!
% Total Variation

alpha = 0.125;

N = n*m;
H = sparse(5*N,5*N);
H(1:N,1:N) = 2*(A.'*A);
h = alpha*ones(5*N,1);
h(1:N) = -2*A.'*md(:);
Aeq         = [Lx,-speye([N,N]),speye([N,N]), sparse(N,N), sparse(N,N);
               Ly, sparse(N,N), sparse(N,N),-speye([N,N]),speye([N,N])];
beq         = zeros(N,2);
lb          = [-Inf(N,1);zeros(4*N,1)];
ub          = Inf(5*N,1);
AA = spdiags([zeros(N,1);-ones(4*N,1)],0, 5*N, 5*N);
iniguess    = zeros(5*N,1);
b           = zeros(5*N,1);
MAXITER = 20; 
QPopt   = optimset('quadprog');
QPopt   = optimset(QPopt,'MaxIter', MAXITER,'Algorithm',...
    'interior-point-convex','Display','iter');
[uvv,val,ef,output] = quadprog(H,h,AA,b,Aeq,beq,lb,ub,iniguess,QPopt);

R4 = uvv(1:N);
r4 = reshape(R4,[n,m]);
figure(6)
imshow(r4);
truesize([400 400])
title(['TV, \alpha = ',num2str(alpha)])
set(gca,'FontSize',20)
