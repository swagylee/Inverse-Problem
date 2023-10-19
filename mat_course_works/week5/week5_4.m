p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);

n = 256;
xvec = linspace(0,2*pi,n);
f = sin(xvec);
f = f(:);

A = convmtx(p,n);
nu = (length(p)-1)/2;
A = A((nu+1):(end-nu),:);

m0 = A*f;
m1 = m0 + 0.1/norm(f)*randn(size(m0));

alphas = [0.0001;0.001;0.01;0.1;1;10;100];

L = eye(n);
L = L-[L(:,end),L(:,1:end-1)];

H = zeros(3*n);
H(1:n,1:n) = 2*(A.'*A);

h = alphas(6)*ones(3*n,1);
h(1:n) = -2*A.'*m1;


Aeq         = [L,-eye(n),eye(n)];
beq         = zeros(n,1);
lb          = [repmat(-Inf,n,1);zeros(2*n,1)];
ub          = repmat(Inf,3*n,1);
AA          = -eye(3*n);
AA(1:n,1:n) = zeros(n,n);
iniguess    = zeros(3*n,1);
b           = [repmat(10,n,1);zeros(2*n,1)];
MAXITER = 200; % Maximum numbers of iterations, Matlab's default value is 200
QPopt   = optimset('quadprog');
%QPopt   = optimset(QPopt,'MaxIter', MAXITER,'Algorithm',...
%    'interior-point-convex','Display','off');
QPopt   = optimset(QPopt,'MaxIter', MAXITER,'display','off');
[uvv,val,ef,output] = quadprog(H,h,AA,b,Aeq,beq,lb,ub,iniguess,QPopt);
recn = uvv(1:n);


% create a subplot for the original signal
subplot(2, 1, 1);
plot(f);
title('Original Signal');
xlabel('Index');
ylabel('Value');

% create a subplot for the reconstructed signal
subplot(2, 1, 2);
plot(recn);
title('Reconstructed Signal by TV Method');
xlabel('Index');
ylabel('Value');