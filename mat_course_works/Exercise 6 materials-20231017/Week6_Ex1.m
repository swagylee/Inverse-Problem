%%% Inverse Problems 1 - 2022
%%% Week 6 - Exercise 1

%%% DON'T REMOVE THESE! OTHERWISE YOU RUN THE RISK OF USING VARIABLES WHICH
%%% ARE NOT DEFINED IN THE CODE.

load('W6_Ex1.mat')

%% %%% Task: Load data stored in 'W6_Ex1.mat', then
%%% (optional but recommended) utilize fftshift in b)-d)
%%% a) Plot the time-series data 'd' in figure 1 with appropriate titles etc.
%%% 
figure(1);
plot(d);
title('Original Signal');


%%% b) Compute the FFT of 'd' and plot the real component and imaginary component 
%%% in figure 2 (with appropriate titles etc.).
%%% 
D = fft(d);
D = fftshift(D);

real_D = real(D);
imag_D = imag(D);


figure(2);

subplot(2, 1, 1);
plot(real_D, 'b');
subplot(2, 1, 2);
plot(imag_D, 'r');

title('Real and  imaginary component');



%%% c) Plot the magnitude of the FFT in figure 3.
%%% 
figure(3);
magnitude_D = abs(D);
plot(magnitude_D);

%%% Plot the phase angle of the FFT coefficients. The function 'angle'
%%% might be useful. Also consider only plotting those frequencies which
%%% are "meaningful".

important_indices_20 = find(abs(D) > 20);
phase_20 = angle(D(important_indices_20));
num_points_20 = length(important_indices_20);

important_indices_10 = find(abs(D) > 10);
phase_10 = angle(D(important_indices_10));
num_points_10 = length(important_indices_10);

figure(4);
plot(phase_20);
% Scatterplot is a good choice here but feel free to try different things
% scatter(...)
% These might be helpful to adjust the y ticks.
ylim = [-pi,pi];
yticks(pi*linspace(-1,1,5))
yticklabels({'-\pi', '-\pi/2', 0, '\pi/2', '\pi'})
