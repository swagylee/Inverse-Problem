p = [1;1;1];
p = conv(p,p);
p = p/sum(p);

nu = (length(p)-1)/2;
n = 10;

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


svd(B)

