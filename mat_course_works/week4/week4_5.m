p = [1; 1; 1];
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);


nu = (length(p)-1)/2;
n = 256;

A = convmtx(p,n);
A = A((nu+1):(end-nu),:);


g = @(x) (2<x & x<=5).*(1)+(5<x & x<=8).*(0.5)+...
    (8<x & x<=11).*(0.75)+(11<x & x<=14).*(0.25);
x1 = linspace(0,16,256);
f = g(x1)';

rng(0,'twister');
noise = 1/norm(f)*randn(size(A*f));
m0 = A*f;
m1 = m0 + 0.1*noise;
m2 = m0 + 0.5*noise;
m3 = m0 + 0.8*noise;
m4 = m0 + noise;

alphas = 10.^[-8:0.1:2];

L = eye(size(A));

Lx_norms = zeros(length(alphas),1);
Ax_m_norms = zeros(length(alphas),1);
L2_errors = zeros(length(alphas),1);
for ii = 1:length(alphas)
    m = m4;
    alpha = alphas(ii);
    syst  = [A; sqrt(alpha)*L];
    m_aug = [m(:);zeros(length(m),1)];
    rec = syst\m_aug;
    Lx_norms(ii)   = norm(L*rec);
    Ax_m_norms(ii) = norm(A*rec-m(:));
    L2_errors(ii) = norm(rec(:)-f(:),2)/norm(f(:),2)*100;

    figure(1)
    clf
    hold off
    plot(f,'k')
    hold on
    plot(rec,'r')
    axis([1 128 -1.2 1.2])
    axis square
    box off
    drawnow
end

figure(2)
endind = round(.8*length(Ax_m_norms));
plot(log(Ax_m_norms(1:endind)),log(Lx_norms(1:endind)),'k')

log_lx = log(Lx_norms(1:endind));
