p = [1; 3; 8; 3; 1];
p_normalized = p / sum(p);
len_f = 16;
M = (length(p)-1)/2;
conv_matrix = convmtx(p_normalized, len_f);

conv_matrix = conv_matrix((M+1):(end-M),:);
%conv_matrix = conv_matrix(:,(M+1):(end-M));

B = conv_matrix;
B(1:M, :) = B(1:M, :) + conv_matrix(end-M+1:end, :);
B(end-M+1:end, :) = B(end-M+1:end, :) + conv_matrix(1:M, :);

f = [0;0;0;0;1;1;0;0;0;0;0;0;0;0;1;2];  % Defined vector

conv_result = conv(f, p_normalized, 'same');

% Calculate the convolution using the convolution matrices A and B
Af = conv_matrix * f;
Bf = B * f;

isequal(Af, conv_result)  % Should be true if A*f is equal to conv(f, p, 'same')
isequal(Bf, conv_result)  % Should be true if B*f is equal to conv(f, p, 'same')
isequal(Af, Bf)  % Should be true if A*f is equal to B*f
