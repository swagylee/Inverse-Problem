% Example computations related to one-dimensional deconvolution.
% Demonstration of sparsity-based choice of regularization parameter for
% total variation regularized reconstructions.
%
% This routine needs Optimization Toolbox (because we use quadprog.m).
%
% Routines DC02_discretedata_comp.m and DC03_nocrimedata_comp.m must be
% precomputed before this file.
%
% Samuli Siltanen February 2014

% Noise level
sigma1 = 0.05;

% Collection of regularization parameters
alphavec = 10.^linspace(-6,4,10);
    
% Plot parameters
lwidth = .5;
thickline = 2;
fsize = 20;
msize = 4;

% Load data computed using a very fine grid, avoiding inverse crime
%load DC2_discretedata A x xx n m mn mIC sigma
%%%%%%%%
% Create your signal, matrix and data here
% Define a signal f (and x to calculate values at)
% A is the convolution matrix
% m is the measurement m = A*f
% mn is the noisy measurement
% x is the discretization of the signal

recx = x;
% Choose data
data = mn(:);
n = length(data);


% Construct prior matrix of size (n)x(n). This implements difference
% between consecutive values assuming periodic boundary conditions.
L = eye(n);
L = L-[L(:,end),L(:,1:end-1)];

% Construct input arguments for quadprog.m
H           = zeros(3*n);
H(1:n,1:n)  = 2*A.'*A;
Aeq         = [L,-eye(n),eye(n)];
beq         = zeros(n,1);
lb          = [repmat(-Inf,n,1);zeros(2*n,1)];
ub          = repmat(Inf,3*n,1);
AA          = -eye(3*n);
AA(1:n,1:n) = zeros(n,n);
b           = [repmat(10,n,1);zeros(2*n,1)];
iniguess    = zeros(3*n,1);

% Set maximum numbers of iterations
MAXITER = 800; % Matlab's default value is 200
QPopt   = optimset('quadprog');
QPopt   = optimset(QPopt,'MaxIter', MAXITER);

% Create plot window and plot original signal
figure(1)
clf
plot(xx,DC_target(xx),'k','linewidth',lwidth)
hold on

% Loop over regularization prameters and compute sparsity of
% reconstructions
sparsevec = zeros(size(alphavec));
recomat   = zeros(n,length(alphavec));
for iii = 1:length(alphavec)
    
    alpha = alphavec(iii);
    
    % Compute regularized solution by constrained quadratic programming
    % using alpha as regularization parameter
    f = [-2*A.'*data; repmat(alpha,2*n,1)];
    [uvv,val,ef,output] = quadprog(H,f,AA,b,Aeq,beq,lb,ub,iniguess,QPopt);
    rec = uvv(1:n);
    disp(['Number of iterations: ', num2str(output.iterations)])
    
    % Record the reconstruction
    recomat(:,iii) = rec(:);
    
    % Take a look at the reconstruction
    plot(recx,rec,'k','linewidth',lwidth)
    
    % Axis settings
    axis([0 1 -.5 2])
    set(gca,'xtick',[0 1/2 1],'fontsize',fsize)
    set(gca,'ytick',[0  1 2],'fontsize',fsize)
    set(gca,'xticklabel',{})
    set(gca,'yticklabel',{})
    set(gca,'PlotBoxAspectRatio' , [2 1 1])
    box off
    
    % Compute sparsity of reconstruction
    sparsevec(iii) = sum(double(abs(L*rec(:))>1e-9));

    % Monitor the run
    disp([' '])
    disp(['Done ', num2str(iii), ' out of ', num2str(length(alphavec))])
    disp([' '])
end

% Compute sparsity of original signal
xxx     = linspace(0,1,length(recx));
origsig = DC_target(xxx);
truesparsity = sum(double(abs(L*origsig(:))>1e-9));

% Save results to file
save DC09_sparsity truesparsity recx recomat sparsevec alphavec
