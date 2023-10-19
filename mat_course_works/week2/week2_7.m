p = [1;5;0.5];
p = conv(p,p);
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);

M = (length(p)-1)/2;

A = convmtx(p,1024);

A = A((M+1):(end-M),:);

d = svd(A);

semilogy(d);

largest_singular_value = d(1);
smallest_singular_value = d(end);

norm(A)

cond(A)