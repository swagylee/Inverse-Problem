img = imread('logo.tif');

f = im2double(img);
figure(1);
imshow(f);
n = 107;
m = 122;
px = 3;
py = 5;
nux = (px-1)/2;
nuy = (py-1)/2;
P = fspecial('gaussian',[px,py],2);

A = convmtx2(P,n,m);
Nc = (n+px-1)*(m+py-1);
I = [1:(n+px-1)*nuy,(Nc-nuy*(n+px-1)+1):Nc];
for i = 1:nux
    I = [I, i:(n+px-1):Nc, (n+px-i):(n+px-1):Nc];
end
A(I,:) = [];

f = f(:);
f_blurred_vec = A*f; %calculate the convolution here

f_blurred = reshape(f_blurred_vec,[n,m]);%reshape the result back to original size
figure(2);
imshow(f_blurred);

mu = 0;
delta = 0.01;
md = imnoise(f_blurred, 'gaussian',mu,delta);
figure(3);
imshow(md);

P = fspecial('motion',3,45);

f = im2double(img);
f_motion_blurred = conv2(f,P,'same');

f_noisy = imnoise(f_motion_blurred, 'salt & pepper');
figure(4);
imshow(f_motion_blurred);

