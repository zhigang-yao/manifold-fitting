% xiayq @ 5/23/2022
% xiayq0121@zufe.edu.cn
% refered to Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228                                                                                                                                                                                                                                                                                                                                                                          clear; % clc

% create a log file
logname = 'out/log_ours.txt';
% select parameters
opts.epsilon = 1e-8;
opts.maxiter = 10;
opts.display = 1;
opts.beta = 2;
opts.logname = logname;

paras = [10;10];
fp = fopen(logname,'a');

for k = 1 : size(paras,2)
    K = paras(1,k);
    d = paras(2,k);
    

% parameter setting
for rate = [0.2 0.3 0.4]% 0.1 0.5]

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
test_idx = [50, 300, 740, 850, 970];
n = length(test_idx);
all_idx = ones(1,N);
all_idx(test_idx) = 0;
train_idx = find(all_idx);




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

r = mean(Dist2(K,:));



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




