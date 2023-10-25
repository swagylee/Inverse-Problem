function H = haarmtx(n,levels)
% Creates the matrix corresponding to the Haar wavelet transfrom in 1D
% It works only if n is a power of 2 and levels is less than log_2(n)
% The columns of H contain the wavelet basis vectors, ordered as follows:
% the first 2^(L-levels) columns contain the father wavelets on the coarsest level
% the following 2^(L-levels) columns contain the mother wavelets on the coarsest level
% the following 2^(L-levels+1) columns contain the mother wavelets on the second coarsest level
% ...
% the last 2^(L-1) columns contain the mother wavelets on the finest level

L = round(log(n)/log(2));

% Requirements
if 2^L ~= n
    error('The size of the matrix must be a power of 2')
end

if nargin == 1
    levels = L;
end

if levels > L
    error('Too many levels required')
end

% Create the matrix
H=zeros(n);

% Mother wavelets
for j=1:levels
    s = 2^j; % size of the support of a j-th level wavelet
    nws = n/s; % number of wavelet in this level
    for k=1:nws
        v = zeros(n,1); % create the column associated to \phi_{j,k}
        first = (k-1)*s+1; % first index of positive support
        half  = (k-1/2)*s; % last index of positive support
        last  = k*s; % last index of the (negative) support
        v(first:half) = 1;
        v(half+1:last) = -1;
        col_idx = 2^(L-j)+k; % index of wavelet as a column of H
        H(:,col_idx) = v/norm(v); % normalization
    end
end

% Father wavelets
s = 2^levels; 
nws = n/s; % if levels = L, nws = 1;
for k=1:nws
    v = zeros(n,1);
    first = (k-1)*s+1; % beginning of the support
    last = k*s; % end of the support
    v(first:last) = 1;
    col_idx = k; % stored in the first columns
    H(:,col_idx) = v/norm(v);
end

end
