img = imread('logo.tif');

f = im2double(img);

px = 5;
py  = 5;
nux = (px-1)/2;
nuy = (py-1)/2;
P = fspecial('gaussian',[px,py],2);

measurement = imnoise(conv2(f,P,'same'), 'gaussian',0,0.01);

n = 107;
m = 122;
A = convmtx2(P,n,m);
Nc = (n+px-1)*(m+py-1);
I = [1:(n+px-1)*nuy,(Nc-nuy*(n+px-1)+1):Nc];
for i = 1:nux
    I = [I, i:(n+px-1):Nc, (n+px-i):(n+px-1):Nc];
end

A(I,:) = [];


r_nav = reshape((A\measurement(:)),[n,m]);
r_wnr = deconvwnr(measurement,P,0.01/var(measurement(:)));
r_reg = deconvreg(measurement,P,0.01*n*m);

% create a subplot for the original signal
subplot(2, 1, 1);
imshow(f);
title('Original Signal');

% create a subplot for the reconstructed signal
subplot(2, 1, 2);
imshow(r_reg);
title('Reconstructed Signal');

