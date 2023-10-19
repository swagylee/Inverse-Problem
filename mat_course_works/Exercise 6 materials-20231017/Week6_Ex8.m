img = imread('logo.tif');

f = im2double(img);

px = 5;
py  = 5;
nux = (px-1)/2;
nuy = (py-1)/2;
P = fspecial('gaussian',[px,py],2);

n = 107;
m = 122;

Lx = spdiags([ones(n*m,1) -ones(n*m,1)], [0,-n], n*m, n*m);
Ly = spdiags([ones(n*m,1) -ones(n*m,1)], [0,1], n*m, n*m);
Ly((m:m:n*m-1)+1,(m:m:n*m)) = 0;

A = convmtx2(P,n,m);
Nc = (n+px-1)*(m+py-1);
I = [1:(n+px-1)*nuy,(Nc-nuy*(n+px-1)+1):Nc];
for i = 1:nux
    I = [I, i:(n+px-1):Nc, (n+px-i):(n+px-1):Nc];
end
A(I,:) = [];


alphas = [0.01 0.05 0.1 0.25 0.5];

mu = 0;
delta = 0.01;
rng(0,'twister')
md = imnoise(conv2(f,P,'same'), 'gaussian',mu,delta);

for ii=1:length(alphas)
    alpha = alphas(ii);
    N = 13054;
    
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
    'interior-point-convex','Display','off');
    
    [uvv,val,ef,output] = quadprog(H,h,AA,b,Aeq,beq,lb,ub,iniguess,QPopt);
    
    rec = uvv(1:N);

    err = norm(rec(:)-f(:))/norm(f(:))

end
