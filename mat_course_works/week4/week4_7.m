p = [1;1;1];
for iii = 1:4
    p = conv(p,p);
end
p = p/sum(p);


nu = (length(p)-1)/2;
n = 1024;

A = convmtx(p,n);
A = A((nu+1):(end-nu),:);

g = @(x) (5*pi/6<x & x<=2*pi).*(.5*x-3.5)+...
    (2*pi<x & x<=10).*(cos(3.*x))+(10<x & x<=14).*(0.25);
x1 = linspace(0,16,n);
f = g(x1)';

m = A*f;
delta = 1e-1;
rng(8,'twister');
noise = randn(n,1);
noise = noise/norm(noise);
md = m + delta*noise;
norm(delta*noise)
ALPHA = logspace(-9,9,200);

L = eye(size(A));
L2_errors = zeros(length(ALPHA),1);

for ii = 1:length(ALPHA)
    m = md;
    alpha = ALPHA(ii);
    syst  = [A; sqrt(alpha)*L];
    m_aug = [m(:);zeros(length(m),1)];
    rec = syst\m_aug;
    L2_errors(ii) = norm(rec(:)-f(:),2)/norm(f(:),2)*100;
end
residual = zeros(length(ALPHA),1);
index=0
for ii = 1:length(ALPHA)
    m = md;
    alpha = ALPHA(ii);
    syst  = [A; sqrt(alpha)*L];
    m_aug = [m(:);zeros(length(m),1)];
    rec = syst\m_aug;
    residual(ii) = norm(A*rec-md);
    if residual(ii) == norm(delta*noise)
        index = ii
    end
end

Lx_norms = zeros(length(ALPHA),1);
Ax_m_norms = zeros(length(ALPHA),1);

for ii = 1:length(ALPHA)
    m = md;
    alpha = ALPHA(ii);
    syst  = [A; sqrt(alpha)*L];
    m_aug = [m(:);zeros(length(m),1)];
    rec = syst\m_aug;
    Lx_norms(ii)   = norm(L*rec);
    Ax_m_norms(ii) = norm(A*rec-m(:));
end

figure(2)
endind = round(.8*length(Ax_m_norms));
plot(log(Ax_m_norms(1:endind)),log(Lx_norms(1:endind)),'k')

log_lx = log(Lx_norms(1:endind));
m = md;
m_aug = [m(:);zeros(length(m),1)];

alpha = ALPHA(67);
alpha_1 = ALPHA(71);
alpha_2 = ALPHA(26);
syst  = [A; sqrt(alpha)*L];
rec = syst\m_aug;

syst1  = [A; sqrt(alpha_1)*L];
rec1 = syst1\m_aug;

syst2  = [A; sqrt(alpha_2)*L];
rec2 = syst2\m_aug;
% create a subplot for the original signal
subplot(4, 1, 1);
plot(f);
title('Original Signal');
xlabel('Index');
ylabel('Value');

% create a subplot for the reconstructed signal
subplot(4, 1, 2);
plot(rec);
title('Reconstructed Signal non morokov');
xlabel('Index');
ylabel('Value');

% create a subplot for the reconstructed signal
subplot(4, 1, 3);
plot(rec1);
title('Reconstructed Signal  morokov');
xlabel('Index');
ylabel('Value');

% create a subplot for the reconstructed signal
subplot(4, 1, 4);
plot(rec2);
title('Reconstructed Signal L-curve');
xlabel('Index');
ylabel('Value');
