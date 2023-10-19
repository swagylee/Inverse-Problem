img = imread('aurinko_paint2.png');
im = img(:,:,1);
f = im2double(im);

px = 5;
py  = 5;
nux = (px-1)/2;
nuy = (py-1)/2;
P = fspecial('gaussian',[px,py],1);

mu = 0;
delta = 0.001;
rng(0,'twister')
md = imnoise(conv2(f,P,'same'), 'gaussian',mu,delta);

n = 128;
m = 128;

A = convmtx2(P,n,m);
Nc = (n+px-1)*(m+py-1);
I = [1:(n+px-1)*nuy,(Nc-nuy*(n+px-1)+1):Nc];
for i = 1:nux
    I = [I, i:(n+px-1):Nc, (n+px-i):(n+px-1):Nc];
end
A(I,:) = [];

Lx = spdiags([ones(n*m,1) -ones(n*m,1)], [0,-n], n*m, n*m);
Ly = spdiags([ones(n*m,1) -ones(n*m,1)], [0,1], n*m, n*m);
Ly((m:m:n*m-1)+1,(m:m:n*m)) = 0;

F0 = [zeros(floor(n/2)-1,m);0.5*ones(1,m); ones(ceil(n/2),m)];

alphas = [0.01 0.05 0.1 0.5 1 5];


for ii=1:length(alphas)

    alpha = alphas(ii);
    R3 = (A'*A + alpha*(Lx'*Lx) + alpha*(Ly'*Ly))\(A'*md(:)+(alpha*(Lx'*Lx) + alpha*(Ly'*Ly))*F0(:));
    err = norm(R3(:)-f(:))/norm(f(:))
    r_genTik = reshape(R3,[n,m]);
    figure(ii)
    imshow(r_genTik)

end

