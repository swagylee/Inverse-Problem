% This script computes the L1 and L2 norm of the function h(x) defined in
% the example 1 of IP_9. This gives an idea of the penalization done by the
% two different norms.

% Rashmi Murthy October 2019.


% Define the interval on which the piecewise constant function is defined.

x = linspace(0,1,101);

% Define the piecewise constant function

h = zeros(1,length(x));


for ii = 1:length(x)
   if (x(ii) >0) && (x(ii) <= 0.5) 
       f(ii) = 1;
   elseif (x(ii)>0.5)&&(x(ii)<1)
       f(ii) = -1;
   end
    
end


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

disp(L_1)
disp(L_2)


