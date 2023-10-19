p = [1;1;1];
for iii = 1:4
    p = conv(p,p);
end
p = p/sum(p);

nu = (length(p)-1)/2;
n = 128

A = convmtx(p,n);
A = A((nu+1):(end-nu),:);

p_temp = zeros(1,n); 
B = zeros(n); 
p_temp(1:nu+1) = p(nu+1:end);
p_temp(end-(nu-1):end) = p(1:nu);
for jj = 1:n
    B(jj,:) = p_temp(1:n);
    p_temp = circshift(p_temp,1); 
end

figure(1)
subplot(121)
spy(A)
subplot(122)
spy(B)

svd(A)
svd(B)