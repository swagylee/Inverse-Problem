% To find the solution of the system of equations in example 1


%  Define The matrix A

A = [ -1 1; -2 1];

% Define the vector m

m = [0;3];


% the system given to us is Af = m. We need to find f

% By Naive inversion

f = inv(A)*m;

% The components of the vector f are

x = f(1)
y = f(2)

% Define the two linear equations of example 1

xplot =  linspace(-5,5,11);
yplot_1 = xplot;
yplot_2 = 2*xplot+3;

% Plot the two lines

figure(1)
plot(xplot,yplot_1,'b','LineWidth',2)
hold on 
plot(xplot,yplot_2,'r','LineWidth',2)
plot(x,y,'k.')

