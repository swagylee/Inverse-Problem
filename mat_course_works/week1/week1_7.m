first_group = [10, 30, 50, 40, 70, 60, 70, 50, 40, 30];  % Ages of the first group
second_group = [8,32,72,31,43,62,68,42,45,33];
% Define a moving average filter
filter = ones(1, 3) / 3; 
result1 = zeros(1,length(first_group));
result2 = zeros(1,length(second_group));
for i=1:length(first_group)
    if i == 1
        j = 10;
    else
        j = i - 1;
    end
    if i == 10
        h = 1;
    else
        h = i + 1;
    end
    sum_up1 = first_group(i) + first_group(j) + first_group(h);
    sum_up2 = second_group(i) + second_group(j) + second_group(h);
    result1(i)=sum_up1/3;
    result2(i)=sum_up2/3;
end

supremum_norm = max(abs(result1 - result2));

% Print the supremum norm with four decimal places
fprintf('Supremum norm: %.4f\n', supremum_norm)
    
