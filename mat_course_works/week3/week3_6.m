p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);



N = 128;


n = N;
M = (length(p)-1)/2;

A = convmtx(p,n);
A = A((M+1):(end-M),:);


g = @(x) (2<x & x<=5).*(1)+(5<x & x<=8).*(0.5)+...
    (8<x & x<=11).*(0.75)+(11<x & x<=14).*(0.25);
x1 = linspace(0,16,N);
f = g(x1)';
plot(f)
[U,D,V] = svd(A);
d = diag(D);

r_alphas = [10 30 50 72 90];

m0 = A*f;
m1 = A*f+0.01/norm(f)*randn(size(A*f));
m2 = A*f+0.05/norm(f)*randn(size(A*f));
m3 = A*f+0.1/norm(f)*randn(size(A*f));
m4 = A*f+0.5/norm(f)*randn(size(A*f));   

w = zeros(1000,1);

for hh = 1:1000
    m0 = A*f;
    m1 = A*f+0.01/norm(f)*randn(size(A*f));
    m2 = A*f+0.05/norm(f)*randn(size(A*f));
    m3 = A*f+0.1/norm(f)*randn(size(A*f));
    m4 = A*f+0.5/norm(f)*randn(size(A*f)); 
    max_error = inf;
    worst_r_alpha = 0;
    m = m4;
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
