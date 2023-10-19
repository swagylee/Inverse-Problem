%%% Inverse Problems 1 - 2022
%%% Week 6 - Exercise 4

%%% DON'T REMOVE THESE! OTHERWISE YOU RUN THE RISK OF USING VARIABLES WHICH
%%% ARE NOT DEFINED IN THE CODE.
clear all
close all

%% Trying out naive FFT deconvolution
%%% a) Create a sawtooth function of length L = 128 and visualize it.

L = 128;

fs=32; %sampling frequency in Hz

Ts=1/fs; %sampling period in s

f=1; %fundamental frequency of the signal in Hz

T=1/f; % perid of the signal in s

nPeriods=4;

Npoints=nPeriods*T/Ts;

time_vector=linspace(0,nPeriods*T-Ts,Npoints); %time vector in seconds

x = sawtooth(2*pi*f*time_vector);

figure();

plot(time_vector,x)

grid on



%%% Note: Matlab has built-in sawtooth function with period of 2pi.  

%%% b) Create Gaussian convolution kernel (in the spatial domain) and use
%%% it to smooth the original signal (keep it length L). Plot the smoothed
%%% version as well.

psf = gausswin(13);

psf = psf/sum(psf);

m = conv(x,psf,'same');

figure(2);

plot(time_vector,m)

grid on

%%% Note: Normalizing the gaussian makes life easier

%%% c) Try deconvolution by computing the (length L) FFT of the PSF and
%%% then pointwise divide the FFT of the smooth signal by the FFT of the PSF,
%%% IFFT and see if you can get the original function (or something similar) back.

W = fft(psf,L);
W_shifted = fftshift(W);
figure(3);

plot(W_shifted)

grid on

% Compute the FFT of m
M = fft(m, L);
M_shifted = fftshift(M);

% Compute the maximum and minimum values of |M|
max_M = max(abs(M_shifted));
min_M = min(abs(M_shifted));

% Display the maximum and minimum values of |M|
fprintf('Max of |M|: %f\n', max_M);
fprintf('Min of |M|: %f\n', min_M);

% Plot the |M|
figure(4);
plot(abs(M_shifted));
title('FFT of the Smoothed Signal');
grid on;


D = M_shifted ./ W_shifted;

% Return to the spatial domain using IFFT
d = ifft(D);

% Plot the result
figure(5);
plot(real(d));
title('Result of Naive Deconvolution');
grid on;





