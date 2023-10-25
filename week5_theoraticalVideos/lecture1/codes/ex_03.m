% This script computes the L1 and L2 norm of the function f(x) defined in
% the example 2 of IP_9. This gives an idea of the penalization done by the
% two different norms.

% Rashmi Murthy October 2019.


% Define the interval on which the piecewise constant function is defined.

x = linspace(0,1,101);

% Define the piecewise constant function

f = zeros(1,length(x));


   f = sin(2*pi*x);



N = length(x)-1;

x = x(1:N);
f = f(1:N);

figure(1)
clf
plot(x,f,'k')



% Define the Differentiation operator

L = -eye(N)+diag(ones(1,N-1),1);
L = L(1:N-1,:);

% Compute Lh

Lh = 10.*(L*f(:));


% Compute the one norm and two norm0

L_1 = norm(Lh,1);
L_2 = (norm(Lh,2))^2;

display(L_1)
disp(L_2)


