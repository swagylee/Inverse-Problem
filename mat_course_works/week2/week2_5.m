p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);

M = (length(p)-1)/2;

A = convmtx(p,1024);

A = A((M+1):(end-M),:);

d = svd(A);

norm(inv_A);

det(A)