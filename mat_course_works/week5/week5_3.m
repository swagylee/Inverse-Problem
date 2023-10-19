p = [1/16,3/16,1/2,3/16,1/16].';

f = [0,0,1,1,1,0.5,0.5,0.5,0.75,0.75,0.75,0.25,0.25,0.25,0,0].';
n = length(f);

A = convmtx(p,n);

nu = (length(p)-1)/2;

A = A((nu+1):(end-nu),:);


m0 = A*f; % noisless observations
m1 = m0 + 0.1/norm(f)*randn(size(m0));
rng(0,'twister'); % noisy observations


H = zeros(3*n); % y in slides [f,v_+,v_-]
H(1:n,1:n) = 2*(A.'*A); % H in slides

alphas = [0.0001,0.001,0.01,0.1,1,10,100];

L = eye(n);
L = L-[L(:,end),L(:,1:end-1)];

L2errors = zeros(1,length(alphas));

for ii = 1:length(alphas)

    h = alphas(ii)*ones(3*n,1);
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
    
    L2errors(ii) = norm(f-recn)/norm(f);

end



