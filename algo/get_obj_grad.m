function [flag, f, G] = get_obj_grad(x, P, Ui, d, r, Ftype, beta)
% calculate the objective value f, that is, ||Pi_x F(x)||_2^2
% and the gradient G of ||Pi_x F(x)||_2^2
%
% input variables
% P : points for epsilon-net
% data : noisy samples
% x : point to be projected to the manifold
% d : dimension of intrinsic manifold
% r : bandwidth parameter
% Ftype : Fi is selected as Fi = x-pi ('xy') or Fi = Pi(x-pi) ('cf')
%
% output variables
% flag : whether x has enough neighbors or not
% f : objective value
% G : 0.5 times gradient
%
% variables in code
% D : dimension of the ambient space
%
% xiayq @ 8/14/2019
% xiayq@zju.edu.cn
% refered to Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228



D = size(P, 1);

% calculate x - pi
dxp = -bsxfun(@minus, P, x);
dist2 = sum(dxp.^2);

indexing = find(dist2<(r*r));
num_idx = numel(indexing);

flag = num_idx > 0;
if flag == 0
    f = -1; G = zeros(D,1);
    return;
end

%P = P(:,indexing);
dxp = dxp(:,indexing);
dist2 = dist2(:,indexing);
Ui = Ui(:,:,indexing);

% calculate tilde_alpha_i, alpha_i and alpha
tilde_alpha_i = zeros(num_idx,1);
for i = 1:num_idx
    tilde_alpha_i(i) = (1- dist2(i)/(r*r)).^beta;
end
alpha = sum(tilde_alpha_i);
alpha_i = tilde_alpha_i/alpha;

% % calculate Pi
% Pi = zeros(D,D,num_idx);
% for i = 1:num_idx
%     Pi(:,:,i) = getPi(P(:,i),data,d,r,dr_type);
% end

% calculate Ax and Fi = alpha_i Pi(x-pi)
Fi = zeros(D,num_idx);
Ax = zeros(D);
for i = 1:num_idx
    Ax = Ax + alpha_i(i)*(eye(D) - Ui(:,:,i)*Ui(:,:,i)');
    switch Ftype
        case 'cf'
            Pi = eye(D)- Ui(:,:,i)*Ui(:,:,i)';
            Fi(:,i) = Pi*dxp(:,i)*alpha_i(i);
            %Fi(:,i) = (dxp(:,i) - Ui(:,:,i)*(Ui(:,:,i)'*dxp(:,i)))*alpha_i(i);
        case 'xy'
            Fi(:,i) = dxp(:,i)*alpha_i(i);
    end
end
Ax = (Ax+Ax')/2;

% calculate Fx
Fx = sum(Fi,2);

% calculate V and Pix*b
[V, lambda] = eig(Ax);
lambda = diag(lambda);
[lambda, idx] = sort(lambda,'descend');
V = V(:,idx);

% Ps = zeros(D,D,D);
% for i = 1 : D
%     Ps(:,:,i) = V(:,i)*V(:,i)';
% end
% Pix = sum(Ps(:,:,1:D-d), 3);

% calulate the objective value
%b = Pix*Fx;
b = Fx - V(:,D-d+1:end)*(V(:,D-d+1:end)'*Fx);
f = sum(b.^2);
if nargout < 3; return; end

% calculate B and T
B = b * Fx';

T = zeros(D);
v = V(:, 1:D-d);
for j = D-d+1:D
    vj = V(:,j);
    dlambda = 1./(lambda(j)-lambda(1:D-d));
    temp = bsxfun(@times, v, dlambda')*v';
    T = T + temp*(B*vj)*vj'+vj*(vj'*B)*temp;
end


% for i = 1 : D-d
%     for j = D-d+1 : D
%         vi = V(:,i);
%         vj = V(:,j);
%         T = T + 1./(lambda(j)-lambda(i))*...
%             (vi*(vi'*B*vj)*vj'+vj*(vj'*B*vi)*vi');
%     end
% end

% calculate dalpha_i
bar_alphai = tilde_alpha_i.^((beta-1)/beta);
dalpha_i = bsxfun(@times, dxp, bar_alphai');
dalpha_i = dalpha_i - bsxfun(@times, sum(dalpha_i,2), alpha_i');
dalpha_i = -dalpha_i*(2*beta/(r*r*alpha));


% calculate <T, Pi>, <b, Fi>
PiT = zeros(1, num_idx);
for i = 1 : num_idx
    %Pi = eye(D) - Ui(:,:,i)*Ui(:,:,i)';
    %PiT(i) = sum(Pi(:).*T(:));
    ui = Ui(:,:,i);
    Tui = T*ui;
    PiT(i) = trace(T) - sum(ui(:).*Tui(:));
end
Fib = b'*Fi;

switch Ftype
    case 'cf'
        G = sum(bsxfun(@times, dalpha_i, Fib+PiT),2)+Ax*b;
    case 'xy'
        G = sum(bsxfun(@times, dalpha_i, Fib+PiT),2)+b;
end

end
% function Axv = Ax_Ui(Ui,alpha_i,v)
%     [D,~,n] = size(Ui);
%     Axv = zeros(D,1);
%     for i = 1 : n
%         ui = Ui(:,:,i);
%         Axv = Axv + (v - ui*(ui'*v))*alpha_i(i);
%     end
%         
% end

