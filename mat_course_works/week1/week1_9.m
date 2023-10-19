p = [1;1;1];
p = conv(p,p);
p = conv(p,p);
p = p/sum(p);
M = (length(p)-1)/2;
f = zeros(128,1);
f(end/4:end*3/4) = 1;

A = convmtx(p,length(f));
A = A((M+1):(end-M),:);

m0 = A*f;
rng(3,'twister');
noise = 1/norm(f)*randn(size(A*f));
m1 = A*f+0.000001*noise;
rng(4,'twister');
noise = 1/norm(f)*randn(size(A*f));
m2 = A*f+0.001*noise;
rng(7,'twister');
noise = 1/norm(f)*randn(size(A*f));
m3 = A*f+0.01*noise;
rng(9,'twister');
noise = 1/norm(f)*randn(size(A*f));
m4 = A*f+0.1*noise;


% Plot the measurements
figure
plot(m0, 'LineWidth', 2)
hold on
plot(m1, 'LineWidth', 2)
plot(m2, 'LineWidth', 2)
plot(m3, 'LineWidth', 2)
plot(m4, 'LineWidth', 2)
hold off
legend('m0', 'm1', 'm2', 'm3', 'm4')
title('Measurements with Different Levels of Noise')
xlabel('Index')
ylabel('Measurement Value')



% Perform the deconvolution
f0 = inv(A)*m0;
f1 = inv(A)*m1;
f2 = inv(A)*m2;
f3 = inv(A)*m3;
f4 = inv(A)*m4;

% Plot the deconvolved signals
figure
plot(f0, 'LineWidth', 2)
hold on
plot(f1, 'LineWidth', 2)
plot(f2, 'LineWidth', 2)
plot(f3, 'LineWidth', 2)
plot(f4, 'LineWidth', 2)
hold off
legend('f0', 'f1', 'f2', 'f3', 'f4')
title('Deconvolved Signals')
xlabel('Index')
ylabel('Signal Value')


% Calculate the relative errors
relative_error0 = norm(f - f0) / norm(f) * 100
relative_error1 = norm(f - f1) / norm(f) * 100
relative_error2 = norm(f - f2) / norm(f) * 100
relative_error3 = norm(f - f3) / norm(f) * 100
relative_error4 = norm(f - f4) / norm(f) * 100

cond(A)