function [flag,f, G, H] = noiseless_obj_grad_Pi(P, Pi, x, r)
% calculate the objective value f, that is, ||Pi_x F(x)||_2^2
% and the gradient G of ||Pi_x F(x)||_2^2
%
% input variables
% P : points for epsilon-net
% Pi : Projection at each pi, a D*D*n tensor
% x : point to be projected to the manifold
% d : dimension of intrinsic manifold
% r : bandwidth parameter
%
% output variables
% flag : whether x has enough neighbors or not
% f : objective value
% G : 0.5 times gradient
%
% variables in code
% D : dimension of the ambient space
% n : number of points in {pi}
%
% xiayq @ 8/18/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228

[D, n] = size(P);

% calculate the difference bwtween x and each plane say bi=Pi^i(x-pi) 
% and the squared distance fi = ||bi||_2^2, di=sqrt(fi)/2r
bi = zeros(D, n);
dxp = -bsxfun(@minus, P, x);
for i = 1 : n
    bi(:,i) = Pi(:,:,i)*(x-P(:,i));
end
fi = sum(bi.^2);
dists = sum(dxp.^2);
di = sqrt(dists)/(2*r);
indexing = find(di <= 1);
num_idx = numel(indexing);

flag = num_idx > 0;
if flag == 0
    f = -1; 
    G = zeros(D,1);
    H = zeros(D,D);
    return;
end

bi = bi(:,indexing);
fi = fi(:,indexing);
di = di(:,indexing);
dists = dists(:,indexing);
Pi = Pi(:,:,indexing);
dxp = dxp(:,indexing);

% calculate tilde_alpha_i, alpha_i, and alpha
tilde_alpha_i = fun_theta(di);
alpha = max(sum(tilde_alpha_i),eps);
alpha_i = tilde_alpha_i/alpha;

% calculate the objective value
f = sum(fi.*alpha_i);
if nargout < 3; return; end

% calculate dalpha_i
d_tilde_alpha_i = (1/(2*r))./ sqrt(dists) .* fun_dtheta(di);
d_tilde_alpha_i = bsxfun(@times, dxp, d_tilde_alpha_i);
da = sum(d_tilde_alpha_i, 2);
dalpha_i = 1/alpha*(d_tilde_alpha_i - da*alpha_i);
G = sum(bsxfun(@times, dalpha_i, fi),2)+2*sum(bsxfun(@times,bi,alpha_i),2);
if nargout < 4; return; end

% calculate Hessian of tilde_alpha_i
Hi = zeros(D, D, num_idx);
Ax = zeros(D,D);
for i = 1 : num_idx
    b = dxp(:,i);
    a1 = -dists(i)^(-1.5)*fun_dtheta(di(i));
    a2 = 1/(2*r)/dists(i)*fun_ddtheta(di(i));
    a3 = dists(i)^(-0.5)*fun_dtheta(di(i));
    Hi(:,:,i) = ((a1+a2)*(b*b')+a3)/(2*r);
    Ax = Ax + alpha_i(i)*Pi(:,:,i);
end
sum_Hi = sum(Hi,3);

H = zeros(D,D);
for i = 1 : num_idx
    a = dalpha_i(:,i); b = bi(:,i);
    H = H + 2*(a*b'+b*a')-fi(i)/alpha*(da*a'+a*da'-Hi(:,:,i)+alpha_i(i)*sum_Hi);
end
H = H + 2*Ax;
%disp(H);


