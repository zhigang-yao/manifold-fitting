function Ui = getUi(P, data, d, r, dr_type)
% get Pi = UiUi' for each pi
% each column of P is a point in X1, the epsilon-net
% data are the sampled data, denoted by X0 in the paper.
%
% xiayq @ 8/14/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228


[n,N] = size(P);
Ui = zeros(n,d,N);

% rescale data with 1/r;
P = P/r; data = data/r;

for i = 1 : N
    p = P(:,i);
    dis2 = sum(bsxfun(@minus, data, p).^2);
    [dis2,idx] = sort(dis2,'ascend');
    
    id = find(dis2 > 1, 1);
    
    if dis2(1) < eps % delete the point p itself 
        %fprintf('delete p%d from its neighbor\n',i);
        P0 = data(:,idx(2:max(id-1,d+1)));
    else
        P0 = data(:, idx(1:max(id-1,d)));
    end
    
    switch dr_type
        case 'cf1'
            P0 = selectNeighbor(P0,p,d);
            U = orth(bsxfun(@minus,P0,p));
        case 'pca'
            [U,~,~] = svd(bsxfun(@minus, P0, p));
            U = U(:,1:d);
        case 'cf2'
            P0 = selectNeighbor(bsxfun(@minus,P0,p),zeros(size(p)),d);
            U = orth(P0);
    end
    Ui(:,:,i) = U;
end

end

function P0 = selectNeighbor(data, p, d)

[~,N] = size(data);
if N <= d
    P0 = data;
else
    idx = zeros(1,d);
    delta = abs(1-sqrt(sum(bsxfun(@minus,data,p).^2)));
    M = zeros(d,N);
    M(1,:) = delta;
    for i = 1 : d
        [~, curID] = min(max(M,[],1));
        idx(i) = curID;
        M(i,curID) = Inf;
        x = data(:,curID); x = x./norm(x,2);
        if i < d
            M(i+1,:) = abs(sum(x.*data));
        end
    end
    P0 = data(:,idx);
end
end