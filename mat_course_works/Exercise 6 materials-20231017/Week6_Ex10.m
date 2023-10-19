img = imread('sun.png');

im = img(:,:,3);

f = im2double(im); 

n = 242;
m = 242;

px = 5;
py  = 5;
nux = (px-1)/2;
nuy = (py-1)/2;
P = fspecial('gaussian',[px,py],1);

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

alphas = [0.1 0.5 1];
F0 = [zeros(floor(2*n/3),m); ones(ceil(n/3),m)];

for ii=1:length(alphas)
    alpha = alphas(ii);

    % Vectorize the prior knowledge.
    F0_vec = F0(:);
    
    R3 = (A'*A + alpha*(Lx'*Lx) + alpha*(Ly'*Ly))\(A'*f(:));
    % Compute the right-hand side of the equation with the prior knowledge.
    b = A' * f(:) + alpha * (Lx' * Lx + Ly' * Ly) * F0_vec;

    % Solve for R3.
    R3_1 = (A' * A + alpha * (Lx' * Lx + Ly' * Ly)) \ b;

    r_genTik = reshape(R3,[n,m]);
    r_genTik_withf = reshape(R3_1,[n,m]);
    figure(ii)
    title("No prior")
    subplot(2,1,1)
    imshow(r_genTik)
    
    subplot(2,1,2)
    title("Prior")
    imshow(r_genTik_withf)

end
figure(4)

imshow(f)
