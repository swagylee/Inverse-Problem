p = [1;1;1];
p = p/sum(p);
n = 5;
A = myOwnConvmtx(p,n);
disp(A)

function A = myOwnConvmtx(p,n)
    A = convmtx(p,n);
    M = (length(p)-1)/2;
    A = A((M+1):(end-M),:);
end