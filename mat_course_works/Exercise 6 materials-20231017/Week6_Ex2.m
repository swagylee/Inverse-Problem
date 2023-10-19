%%% Inverse Problems 1 - 2022
%%% Week 6 - Exercise 2

%%% DON'T REMOVE THESE! OTHERWISE YOU RUN THE RISK OF USING VARIABLES WHICH
%%% ARE NOT DEFINED IN THE CODE.
clear all;
close all;


% time span
t0=0;
tf=2; %in seconds

% Here we have opted to only plot the fftshifted versions
shiftFFT=true;

%% a) Create your custom fft function in a separate file and save it 
%%    in the same directory of this file


%% b) f(t)=cos(2*pi*f*t)
close all;

Fs=200; % sampling frequency in Hz
f = 10; % cosine frequency in Hz
A=1.0; % amplitude
% continue....

t = t0:1/Fs:tf;
x = A*cos(2*pi*f*t);
[X_shifted, f_shifted] = my_fft(x, Fs, true); % FFT with shift
[X, f] = my_fft(x, Fs, false); % FFT without shift

figure(1);


subplot(3, 1, 1);
plot(t, x);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Signal');

subplot(3, 1, 2);
plot(f, abs(X));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT without Shift');

subplot(3, 1, 3);
plot(f_shifted, abs(X_shifted));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT with Shift');


%% c) f(t)=cos(2*pi*f1*t)+0.5*sin(2*pi*f2*t)

Fs=200; %in Hz
f1 =10; % cosine frequency in Hz
f2 =40; % sine frequency in Hz
A1=1.0; % amplitude 1
A2=0.5; % amplitude 2

% continue....

t = t0:1/Fs:tf;

x_1 = zeros(size(t));
x_1(t > 0.5 & t < 1) = 2;

[X_shifted_1, f_shifted_1] = my_fft(x_1, Fs, true); % FFT with shift
[X_1, f_1] = my_fft(x_1, Fs, false); % FFT without shift

figure(2);


subplot(3, 1, 1);
plot(t, x_1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Signal');

subplot(3, 1, 2);
plot(f_1, abs(X_1));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT without Shift');

subplot(3, 1, 3);
plot(f_shifted_1, abs(X_shifted_1));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT with Shift');


%% d)square pulse between t1 and t2


Fs=200; % in Hz
t1=0.5; % in s
t2=1.0; % in s
A=2;    % amplitude
% continue....

