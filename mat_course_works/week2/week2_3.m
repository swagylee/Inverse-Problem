A = [-1 -1 2; 0 -1 1; -2 0 2];
b = [2; -5; 4];

kernel_A = null(A);

norm(A*kernel_A);

A_inv = inv(A);

A_pinv = pinv(A);

d = svd(A);
[U,D,V] = svd(A);

V*V';
inv(U)

transpose(U)