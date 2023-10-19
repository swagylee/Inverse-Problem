% Define the vectors
p_tilde = [1,1,1];
f = [0,0,0,0,1,1,0,0,0,0];

% Pad f with periodic boundary conditions
f_padded = [f(end), f, f(1)];

% Preallocate the result vector
result = zeros(1, length(f));

% Compute the convolution
for j = 1:length(f)
    result(j) = sum(p_tilde .* f_padded(j:j+2));
end

% Initialize given vectors
p_tilde = [1,1,1];
f = [2,0,0,0,1,1,0,0,0,0];

% Normalize point spread function
p = p_tilde / sum(p_tilde);

% Create a padded version of f using periodic boundary conditions
f_padded = [f(end), f, f(1)];

% Initialize convolution vector
convolution_result = zeros(1, length(f));

% Perform convolution
for j = 1:length(f)
    for l = -1:1
        convolution_result(j) = convolution_result(j) + p(l+2) * f_padded(j-l+1);
    end
end

% Display the result
disp(convolution_result);