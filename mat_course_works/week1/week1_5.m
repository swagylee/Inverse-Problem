v1 = [1;1;1];
v2 = [0;1;1;1;0];
v3 = [0;1;2;3;4;3;2;1;0];
v4 = [0;1;2;2;3;3;2;1;0];
v5 = [0;1;2;1;0;-1;-2;-1;0];

%figure; stem(v1); title('v1');
%figure; stem(v2); title('v2');
%figure; stem(v3); title('v3');
%figure; stem(v4); title('v4');
%figure; stem(v5); title('v5');

p_times_f = conv(v1, v2);
len_p_times_f = length(p_times_f);

p_times_f_same = conv(v1, v2, 'same');

conv_v2_v1 = conv(v2, v1, 'same');
conv_v1_v2 = conv(v1, v2, 'same');
conv_v2_v3 = conv(v2, v3, 'same');
conv_v3_v1 = conv(v3, v1, 'same');
conv_v1_v4 = conv(v1, v4, 'same');
conv_v4_v1 = conv(v4, v1, 'same');
conv_v5_v2 = conv(v5, v2, 'same');
conv_v5_v3 = conv(v5, v3, 'same');

figure; stem(conv_v2_v1); title('conv_v2_v1');
figure; stem(conv_v1_v2); title('conv_v1_v2');
figure; stem(conv_v2_v3); title('conv_v2_v3');
figure; stem(conv_v3_v1); title('conv_v3_v1');
figure; stem(conv_v1_v4); title('conv_v1_v4');
figure; stem(conv_v4_v1); title('conv_v4_v1');
figure; stem(conv_v5_v2); title('conv_v5_v2');
figure; stem(conv_v5_v3); title('conv_v5_v3');

