% Example 1 To see what l1 and l2 norms do to a function.
% Here the function is piecewise continuous.

% Let's define the values of x over which h is defined.

x = linspace(0,1,11);

% Define the function f

h = zeros(1,length(x));

for ii = 1:length(x)
   if(x(ii)> 1/2) && (x(ii)<=1)
       h(ii) = 2;
   end
    
end

N = length(x)-1;
x = x(1:N);
h = h(1:N);

% figure(1)
% clf
% plot(x,h,'k','LineWidth',1.5)

% Define the operator L

L = -eye(N)+diag(ones(1,N-1),1);
L = L(1:N-1,:);

% Compute Lh

Lh = 10.*(L*h(:));

l_1_norm = norm(Lh,1);
l_2_norm = norm(Lh,2).^2;

disp(l_1_norm)
disp(l_2_norm)
