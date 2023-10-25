clear all
close all
clc
% Example 2 To see what l1 and l2 norms do to a function.
% Here the function is piecewise continuous.

% Let's define the values of x over which h is defined.

x = linspace(0,1,11);

% Define the function f

f = zeros(1,length(x));

f = 2*x;

N = length(x)-1;
x = x(1:N);
f = f(1:N);

 figure(1)
 clf
 plot(x,f,'k','LineWidth',1.5)
 
 

% Define the operator L

L = -eye(N)+diag(ones(1,N-1),1);
L = L(1:N-1,:);

% Compute Lh

Lf = 10.*(L*f(:));

l_1_norm = norm(Lf,1);
l_2_norm = norm(Lf,2).^2;

disp(l_1_norm)
disp(l_2_norm)
