% Sparsity regularization for one-dimensional deconvolution
% using "half-quadratic" optimization.
%
% Arguments:
% A        convolution matrix
% L        prior matrix
% mn       data vector
% alpha    regularization parameter, positive real number
% MAXITER  maximum number of iterations 
%
% This routine needs Matlab's optimization toolbox.
%

function rec = sparsity_recon(A,L,mn,alpha,MAXITER)

data = mn;
n = length(data);

QPopt   = optimset('quadprog');
QPopt   = optimset(QPopt,'MaxIter', MAXITER,'display','off');

% Construct input arguments for quadprog.m
H           = zeros(3*n);
H(1:n,1:n)  = 2*((A.')*A);

%equality constraint
% Lf=(v+) - (v-)   <==>   Lf - (v+) + (v-) =0
Aeq         = [L,-eye(n),eye(n)];
beq         = zeros(n,1);

% box constraint
% -Inf< f < Inf
% 0 <= v+ < Inf
% 0 <= v- < Inf
lb          = [repmat(-Inf,n,1);zeros(2*n,1)];
ub          = repmat(Inf,3*n,1);

% Inequality contraint
% we don't have any constraint of this type.
% We create a obvious condition that is automatically satisfied
% 0*f<10

AA          = -eye(3*n);
AA(1:n,1:n) = zeros(n,n);
b           = [repmat(10,n,1);zeros(2*n,1)];


% measurement component
% h=[-2A^Tm ; \alpha*ones ; \alpha*ones]
f           = [-2*((A.')*data); repmat(alpha,2*n,1)];

%initial guess
iniguess    = zeros(3*n,1);

% Compute minimizer by constrained quadratic programming
% using alpha as regularization parameter
[uvv,val,ef,output] = quadprog(H,f,AA,b,Aeq,beq,lb,ub,iniguess,QPopt);

% extract only the first n components
rec = uvv(1:n);




