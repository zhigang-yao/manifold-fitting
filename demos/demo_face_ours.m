% To denoise facial images and generate Figure 6-8.
%
%
% xiayq @ 8/19/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228

clear; % clc

logname = 'out/log_ours.txt';
% select parameters
opts.epsilon = 1e-8;
opts.maxiter = 10;
opts.display = 1;
opts.beta = 2;
opts.logname = logname;

paras = [200;10];
fp = fopen(logname,'a');

for k = 1 : size(paras,2)
    K = paras(1,k);
    d = paras(2,k);
    

% parameter setting
for rate = [0.2 0.3 0.4]

% load data
try 
    D = load(sprintf('face_data/sub22_part_rate%.1f.mat',rate));
    X = D.X;
    X0 = D.X0;
    h = D.h;
    w = D.w;
    N = size(X,2);
    fprintf('----data loading succeed for rate %.1f----\n', rate);
    fprintf(fp, '----data loading succeed for rate %.1f----\n', rate);
catch
    load face_data/sub22_part.mat X
    [h, w, N] = size(X);
    D = h*w;
    X0 = reshape(X,D,N);
    sigma = mean(X0(:))*rate;
    % add noise
    E = randn(D,N)*sigma;
    X = X0 + E;
    save(sprintf('face_data/sub22_part_rate%.1f.mat',rate),'X','X0','h','w');
end


% partition dataset    
%train_idx = 1:2:1000;
test_idx = [50, 300, 740, 850, 970];
n = length(test_idx);
all_idx = ones(1,N);
all_idx(test_idx) = 0;
train_idx = find(all_idx);
%disp(numel(train_idx));



Dist2 = zeros(N);
for i = 1 : N
    for j = 1 : N
        if i < j
            dist2 = norm(X(:,i)-X(:,j),'fro');
            Dist2(i,j) = dist2;
        end
    end
    Dist2(i,i) = inf;
end
Dist2 = sort(Dist2 + Dist2');

r = ceil(mean(Dist2(K,:)));
%d = 10;


%%
tic;
X_sample = X(:,train_idx);
X_ini = X(:,test_idx);
[Mout, info] = manfit_our(X_sample, d, r, X_ini, opts);
t = toc;

fprintf('----rate=%.1f, r=%d, d=%d, and costs %.1f seconds---- \n',rate, r, d, t);
fprintf(fp,'----rate=%.1f, r=%d, d=%d, and costs %.1f seconds---- \n',rate, r, d, t);

save(sprintf('out/face/new_face_ours_d%d_K%d_rate%.1f.mat', d, K, rate),'Mout','train_idx','test_idx','t');

end

end

fclose(fp);
return;
%% visualize the results
Mout = reshape(Mout, h, w, n);
test_X0 = reshape(X0(:,test_idx), h, w, n);
test_X = reshape(X_ini, h, w, n);
for i = 1 : n
    subplot(1,3,1);
    imshow(test_X0(:,:,i),[]);
    subplot(1,3,2);
    imshow(test_X(:,:,i),[]);
    subplot(1,3,3);
    imshow(Mout(:,:,i),[]);
    title(i)
    pause;
end


