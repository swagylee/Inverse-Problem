n = 1024;
g = @(x) (5*pi/6<x & x<=2*pi).*(.5*x-3.5)+...
    (2*pi<x & x<=10).*(cos(3.*x))+(10<x & x<=14).*(0.25);
x1 = linspace(0,16,n);
f = g(x1)';

p = [1; 1; 1];
p = conv(p,p);
p = conv(p,p);
p = conv(p,p);
p = p/(sum(p));

nu = (length(p)-1)/2;
A = convmtx(p,n);
A = A((nu+1):(end-nu),:);

m = A*f;

rng(5,'twister');
noise = randn(n,1);
noise = noise/norm(noise);

% your original setup code here
...

% Define your L matrices and their labels
L = {eye(size(A)), eye(n)-diag(ones(n-1,1),-1), zeros(n)-diag(ones(n-1,1),-1)+diag(ones(n-1,1),1)};
L_labels = {'L1', 'L2', 'L3'};

delta = [0.001, 0.01, 0.1, 1, 10, 100, 1000]; % noise levels

% Create a new figure
figure;

% Loop over different L matrices
for i = 1:length(L)
    % Loop over different noise levels
    for j = 1:length(delta)
        % Add noise
        md = m + delta(j)*noise;
        
        % Solve the regularized problem
        fd = (A'*A + alpha*L{i})\(A'*md);
        
        % Compute the error
        err = norm(f-fd)/norm(f);
        fprintf('Error for %s and delta = %f: %f\n', L_labels{i}, delta(j), err);
        
        % Create a new subplot for each noise level and L matrix
        subplot(length(L), length(delta), (i-1)*length(delta) + j);
        plot(x1, f, 'b', x1, fd, 'r');
        title(sprintf('Reconstruction with %s and delta = %f', L_labels{i}, delta(j)));
        legend('Original', 'Reconstruction');
        xlabel('x');
        ylabel('f(x)');
    end
end



