% System of eqiuations in example 2

% Define the matrix A

A = [-1 1;-1 1];

% Define the vector m

m = [0;3];

% Lets compute the naive inversion. That is f = A^{-1}*m;

f = inv(A)*m;

% Defining the linear system

xplot = linspace(-5,5,11);
yplot_1 = xplot;   % First linear equation y = x
yplot_2 = xplot + 3; % Second linear equation y = x+3
% % Let's plot the linear system
% figure(1)
% plot(xplot,yplot_1,'b','LineWidth',2)
% hold on
% plot(xplot,yplot_2,'r','LineWidth',2)

% Compute the inversion using psedoinverse.

f1 = pinv(A)*m;
x = f1(1);  % 1st component of the vector f1
y = f1(2);  % Second component of the vector f1

yplot_ls = xplot +1.5; % Line passing through(x,y) with slope 1



% Let's plot the linear system
figure(2)
plot(xplot,yplot_1,'b','LineWidth',2)
hold on
plot(xplot,yplot_2,'r','LineWidth',2)
plot(x,y,'k*')
plot(xplot,yplot_ls,'g--','LineWidth',2)
