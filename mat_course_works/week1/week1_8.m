% Define the point spread function p
p = [1/16, 3/16, 1/2, 3/16, 1/16]';

% Build the convolution matrix A
A_full = convmtx(p, 16);
A = A_full(3:18, :);

% Define the signal f
f = [0, 0, 1, 1, 1, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.25, 0.25, 0.25, 0, 0]';

% Set the seed for the random number generator to ensure reproducibility
rng(0,'twister');

% Compute the noise vector
noise = 1/norm(f) * randn(size(A*f));

% Compute the measurement vectors m0, m1, m2, m3, and m4
m0 = A*f;
m1 = A*f + 0.01*noise;
m2 = A*f + 0.05*noise;
m3 = A*f + 0.1*noise;
m4 = A*f + 0.5*noise;

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

% Perform deconvolution for all measurements
f0 = inv(A)*m0;
f1 = inv(A)*m1;
f2 = inv(A)*m2;
f3 = inv(A)*m3;
f4 = inv(A)*m4;

% Calculate the relative errors
relative_error0 = norm(f - f0) / norm(f) * 100;
relative_error1 = norm(f - f1) / norm(f) * 100;
relative_error2 = norm(f - f2) / norm(f) * 100;
relative_error3 = norm(f - f3) / norm(f) * 100;
relative_error4 = norm(f - f4) / norm(f) * 100;

% Print the relative errors with two decimal places
fprintf('Relative error with m0: %.2f%%\n', relative_error0)
fprintf('Relative error with m1: %.2f%%\n', relative_error1)
fprintf('Relative error with m2: %.2f%%\n', relative_error2)
fprintf('Relative error with m3: %.2f%%\n', relative_error3)
fprintf('Relative error with m4: %.2f%%\n', relative_error4)

cond(A)
