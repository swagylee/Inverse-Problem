%%% Inverse Problems 1 - 2022
%%% Week 6 - Exercise 3

%%% DON'T REMOVE THESE! OTHERWISE YOU RUN THE RISK OF USING VARIABLES WHICH
%%% ARE NOT DEFINED IN THE CODE.
clear all
close all

nPoints=512;
[X,Y] = meshgrid(0:(nPoints-1),0:(nPoints-1));

%%% NOTE: mat2gray is part of Image Processing Toolbox! If you are unable
%%% to install it via the Matlab Add-On Explorer you can easily replace
%%% it with the following commands which produce pretty much the same
%%% output: myImage = double(im)/255;

% Also note that for some of the images (wave X and wave Y) you need to
% either zoom in or make the figure bigger to actually see the small pixel
% values because the outcome is so sparse.

% By default imshow shows those values which are between [0, 1] and any
% values outside that interval are rounded accordingly. Therefore if you
% image contains values from [0, 255], imshow thinks anything beyond 1
% should be white and you get nonsense. You can either scale the image to
% [0,1] yourself OR as a second argument to imshow give new limits. If you
% set the limits to [], imshow automatically chooses the lower bound to the
% smallest input and upper bound to largest which is often the easiest
% solution.
% However there are cases where we for example know that there shouldn't be
% negative values for example. Then it makes sense to manually set the
% lower bound to 0.
% Finally if you still see nothing, the logarithm is very important because
% of the massive difference between the largest and smallest values: a
% linear greascale colormap simply can not show all of the inbetween values
% in a meaningful way because some large values dominate the image.

%% a) Load the images

%% b) Compute the FFT of the following images

% Waves X 
x=imread('wavesXY.png');

X = fft2(x);
X = fftshift(X);
X_mag = abs(X);

% assume the magnitudes of the FFT2 are stored in (double) matrix X_mag
% first, let’s transform the values to logspace to facilitate visualization
X_mag_transform = log(X_mag+1); % we are adding 1 to avoid problems with log(0)
X_mag_image = mat2gray(X_mag_transform); % here we are transforming the matrix to grayscale values
imshow(X_mag_image); % finally show as an image


% Waves Y
%x=imread('wavesY.png');

% Waves XY 
%x=imread('wavesXY.png');

% Waves circle 
%x=imread('wavesCirc.png');

% Rectangle, aligned to X-Y 
%x=imread('rectangle_01.png');

% Rectangle, at 45 degrees
%x=imread('rectangle_02.png');

%% c) Create 2 plots in the same figure
%%% showing the original image (you can use the function imshow), and the magnitude 
%%% of the shifted FFT coefficients.


