p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);

M = (length(p)-1)/2;

A = convmtx(p,128);

A = A((M+1):(end-M),:);

d = svd(A);

norm(inv(A));

det(A)

semilogy(d);
