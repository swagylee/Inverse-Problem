p = [1;1;1];
for iii = 1:4
    p = conv(p,p);
end
p = p/sum(p);

N = 1024;

n = N;
M = (length(p)-1)/2;

A = convmtx(p,n);
A = A((M+1):(end-M),:);


g = @(x) (5*pi/6<x & x<=2*pi).*(.5*x-3.5)+...
    (2*pi<x & x<=10).*(cos(3.*x))+(10<x & x<=14).*(0.25);

x1 = linspace(0,16,N);

f = g(x1)';



[U,D,V] = svd(A);
d = diag(D);

r_alphas = [10 150 300 600 800];

rng(0,'twister');
noise = 1/norm(f)*randn(size(A*f));  
m0 = A*f;
m1 = A*f+0.01*noise;
m2 = A*f+0.05*noise;
m3 = A*f+0.1*noise;
m4 = A*f+0.5*noise; 

for ii = 1:length(r_alphas)
    % Calculate the Dralpha here
    % Reconstruct with TSVD
    % Do reconstructions here
    r_alpha = r_alphas(ii);
    Dralpha = zeros(size(D));
    d_ralpha = d(1:r_alpha);
    Dralpha(1:r_alpha,1:r_alpha) = diag(1./d_ralpha);
    TSVDrec0 = V*Dralpha*U.'*m0;
    TSVDrec1 = V*Dralpha*U.'*m1;
    TSVDrec2 = V*Dralpha*U.'*m2;
    TSVDrec3 = V*Dralpha*U.'*m3;
    TSVDrec4 = V*Dralpha*U.'*m4;

    subplot(5, 5, 1+(ii-1)*5)
    plot(f,'k')
    hold on
    plot(TSVDrec0,'b')
    axis square

    subplot(5, 5, 2+(ii-1)*5)
    plot(f,'k')
    hold on
    plot(TSVDrec1,'b')
    axis square

    subplot(5, 5, 3+(ii-1)*5)
    plot(f,'k')
    hold on
    plot(TSVDrec2,'b')
    axis square

    subplot(5, 5, 4+(ii-1)*5)
    plot(f,'k')
    hold on
    plot(TSVDrec3,'b')
    axis square

    subplot(5, 5, 5+(ii-1)*5)
    plot(f,'k')
    hold on
    plot(TSVDrec4,'b')
    axis square
end

