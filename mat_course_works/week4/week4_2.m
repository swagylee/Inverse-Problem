p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);

N = 128;
g = @(x) (2<x & x<=5).*(1)+(5<x & x<=8).*(0.5)+...
    (8<x & x<=11).*(0.75)+(11<x & x<=14).*(0.25);
x1 = linspace(0,16,N);
f = g(x1)';

A = convmtx(p,N);

nu = (length(p) - 1)/2;

A = A((nu+1):(end-nu),:);

[U,D,V] = svd(A);

svals = diag(D);


alphas = [0.001, 0.01, 0.1, 0.3, 0.7, 0.9, 3];
m0 = A*f;
m1 = A*f+0.01/norm(f)*randn(size(A*f));
m2 = A*f+0.05/norm(f)*randn(size(A*f));
m3 = A*f+0.1/norm(f)*randn(size(A*f));
m4 = A*f+0.5/norm(f)*randn(size(A*f));  


w = zeros(10000,1);

for jj = 1:10000
    max_value = inf;
    m0 = A*f;
    m1 = A*f+0.01/norm(f)*randn(size(A*f));
    m2 = A*f+0.05/norm(f)*randn(size(A*f));
    m3 = A*f+0.1/norm(f)*randn(size(A*f));
    m4 = A*f+0.5/norm(f)*randn(size(A*f));  
    for ii = 1:length(alphas)
        alpha = alphas(ii); %% current alpha in use
        Dplus_alpha = diag(svals./(svals.^2+alpha));
        Tikh_rec = V*Dplus_alpha*U.'*m4;
        L2error_m4 = norm(f-Tikh_rec)/norm(f)*100;
        if max_value >= L2error_m4
            max_value = L2error_m4;
            w(jj) = alpha;
        end
    end
end

sum(w==3)/length(w)
sum(w==0.9)/length(w)
sum(w==0.7)/length(w)
sum(w==0.3)/length(w)
sum(w==0.1)/length(w)
sum(w==0.01)/length(w)
sum(w==0.001)/length(w)
