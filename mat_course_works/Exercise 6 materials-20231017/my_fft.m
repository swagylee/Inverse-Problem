function [X,f] = my_fft(x,Fs,shift)
%MY_FFT 此处显示有关此函数的摘要
%   Input arguments:
%   x (column vector) the signal array 𝑥𝑛
%   Fs  (float) the sampling frequency of the signal 𝐹𝑠
%   shift (boolean) a flag indicating if we want to shift the FFT or not
%   Output arguments:
%   X  (column vector) the FFT of the signal
%   f  (column vector) the frequencies associated with each 𝑋_k
X = fft(x);
N = length(x);
if shift == false
    k_indexes = 0:N-1;
    f = k_indexes * Fs/N;
else
    X = fftshift(X);
    if mod(N,2) == 0
        k_indexes = -N/2:N/2-1; % N even
    else
        k_indexes = -(N-1)/2:(N-1)/2; % N odd
    end
    f = k_indexes * Fs/N;
end

end

