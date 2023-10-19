
p = [0.25; 0.3; 0.5; 1; 0.9; 0.75; 0.5; 0.25; 0.1];
p = p/sum(p);

M = (length(p)-1)/2;

A = convmtx(p,128);

A = A((M+1):(end-M),:);

d = svd(A);

norm_A = norm(A);

cond_A = cond(A);

cond_A

semilogy(d);
