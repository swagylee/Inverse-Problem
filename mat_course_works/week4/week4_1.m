f = [0,0,1,1,1,0.5,0.5,0.5,0.75,0.75,0.75,0.25,0.25,0.25,0,0].';
p = [1/16; 3/16; 1/2; 3/16; 1/16];

nu = (length(p)-1)/2;
n = length(f);

A = convmtx(p,n);
A = A((nu+1):(end-nu),:);

[U,D,V] = svd(A);
svals = diag(D);

alphas = [0.001, 0.01, 0.1, 0.5, 0.9, 1, 10];


m0 = A*f;
m1 = m0+0.01/norm(f)*randn(size(m0));
m2 = m0+0.05/norm(f)*randn(size(m0));
m3 = m0+0.1/norm(f)*randn(size(m0));
m4 = m0+0.5/norm(f)*randn(size(m0));


for ii = 1:length(alphas)
    alpha = alphas(ii); %% current alpha in use
    Dplus_alpha = diag(svals./(svals.^2+alpha));
    Tikh_rec = V*Dplus_alpha*U.'*m3;
    L2error = norm(f-Tikh_rec)/norm(f)*100;
end

w = zeros(1000,1);

for jj = 1:1000
    max_value = -1;
    m0 = A*f;
    m1 = m0+0.01/norm(f)*randn(size(m0));
    m2 = m0+0.05/norm(f)*randn(size(m0));
    m3 = m0+0.1/norm(f)*randn(size(m0));
    m4 = m0+0.5/norm(f)*randn(size(m0));
    for ii = 1:length(alphas)
        alpha = alphas(ii); %% current alpha in use
        Dplus_alpha = diag(svals./(svals.^2+alpha));
        Tikh_rec = V*Dplus_alpha*U.'*m4;
        L2error_m4 = norm(f-Tikh_rec)/norm(f)*100;
        if max_value < L2error_m4
            max_value = L2error_m4;
            w(jj) = alpha;
        end
    end
end

sum(w==10)/length(w)
sum(w==1)/length(w)
sum(w==0.9)/length(w)
sum(w==0.5)/length(w)
sum(w==0.1)/length(w)
sum(w==0.01)/length(w)
sum(w==0.001)/length(w)