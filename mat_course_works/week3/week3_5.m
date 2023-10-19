% Define the point spread function
p = [1/16; 3/16; 1/2; 3/16; 1/16];

f = [0;0;1;1;1;0.5;0.5;0.5;0.75;0.75;0.75;0.25;0.25;0.25;0;0];


n = length(f);
M = (length(p)-1)/2;

A = convmtx(p,n);
A = A((M+1):(end-M),:);


[U,D,V] = svd(A);


d = diag(D);


max_singular_value = max(d);
min_singular_value = min(d);

r_alphas = [1 5 10 12 16];

m1 = A*f+0.01/norm(f)*randn(size(A*f));
m2 = A*f+0.05/norm(f)*randn(size(A*f));
m3 = A*f+0.1/norm(f)*randn(size(A*f));
m4 = A*f+0.5/norm(f)*randn(size(A*f)); 
measurements = {m0, m1, m2, m3, m4};

w = zeros(1000,1);

for hh = 1:1000
    m0 = A*f;
    m1 = A*f+0.01/norm(f)*randn(size(A*f));
    m2 = A*f+0.05/norm(f)*randn(size(A*f));
    m3 = A*f+0.1/norm(f)*randn(size(A*f));
    m4 = A*f+0.5/norm(f)*randn(size(A*f)); 
    max_error = inf;
    worst_r_alpha = 0;
    m = m1;
    for ii = 1:length(r_alphas)
        r_alpha = r_alphas(ii);
        Dralpha = zeros(size(D));
        d_ralpha = d(1:r_alpha);
        Dralpha(1:r_alpha,1:r_alpha) = diag(1./d_ralpha);
        TSVDrec = V*Dralpha*U.'*m;
        error = norm(f - TSVDrec)/norm(f);
        if error <= max_error
            max_error = error;
            worst_r_alpha = r_alpha;
        end
    end
    w(hh) = worst_r_alpha;
end
sum(w==16)

sum(w==12)

sum(w==10)

sum(w==5)

sum(w==1)
