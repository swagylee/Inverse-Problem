n = 500;
x = (1:n)';

f = x+50*ceil(sin(x/5));

p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);

nu = (length(p)-1)/2;
A = convmtx(p,n);
A = A((nu+1):(end-nu),:);

m = A*f;
delta = 1e-2;
noise = randn(n,1);
noise = noise/norm(noise);
md = m + delta*noise;

alpha = 0.1;

f0 = x+25; %% What should you write here?
fdG = (A'*A + alpha*eye(n))\(A'*md + alpha*f0);

% create a new figure
figure;

% create a subplot for the original signal
subplot(2, 1, 1);
plot(f);
title('Original Signal');
xlabel('Index');
ylabel('Value');

% create a subplot for the reconstructed signal
subplot(2, 1, 2);
plot(fdG);
title('Reconstructed Signal');
xlabel('Index');
ylabel('Value');



