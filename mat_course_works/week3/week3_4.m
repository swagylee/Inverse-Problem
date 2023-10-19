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


f = zeros(128,1);
f(round(end/2):end) = 1;


ma = A*f;
mb = B*f;
rng(0,'twister');
noise = 1/norm(f)*randn(size(A*f));
ma1 = ma + 0.01*noise;
mb1 = mb + 0.01*noise;
ma2 = ma + 0.1*noise;
mb2 = mb + 0.1*noise;

tols = [0.01;0.05;0.1;0.15;0.3;0.4;0.9];

signals = {ma1, ma2, mb1, mb2};
matrices = {A, B};

for i = 1:length(matrices)
    for j = 1:length(signals)
        min_error = inf;
        best_tol = 0;
        for k = 1:length(tols)
            frec = pinv(matrices{i}, tols(k)) * signals{j};
            error = norm(f - frec)/norm(f) * 100;
            if error < min_error
                min_error = error;
                best_tol = tols(k);
            end
        end
        min_error

    end
    
end


% Perform SVD
[UA,DA,VA] = svd(A);
[UB,DB,VB] = svd(B);

% Singular vector indices we want to examine
SV_index = [3, 5, 15, 30, 50];

% Plot singular vectors
figure;
for i = 1:length(SV_index)
    % Singular vector index
    index = SV_index(i);

    % Plot singular vectors of A
    subplot(length(SV_index), 2, 2*i-1);
    plot(VA(:,index));
    title(['A: Singular Vector ', num2str(index)]);

    % Plot singular vectors of B
    subplot(length(SV_index), 2, 2*i);
    plot(VB(:,index));
    title(['B: Singular Vector ', num2str(index)]);
end

