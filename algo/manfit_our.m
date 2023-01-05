function [Mout, info] = manfit_our(data, d, r, X, opts)
% input : 
% data - dataset
% d    - dimension of manifold
% r    - radius of neighborhood
% X    - the set of initial points
% opts - optional parameters, can be set as opts=[];
%        maxiter : maxiter of GD
%        epsilon : stopping criterian
%        alpha0 : initial step length
%        display : whether display iteration information or not
%
% output :
% Mout - output manifold
% info - informations of algorithm
%
% xiayq @ 04/21/2022
%       modify comments
%       replace r by 2r in getUi
% xiayq @ 08/14/2019
%       first version
% xiayq0121@zufe.edu.cn
% refered to Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228


t1 = clock;
Ui = getUi(data, data, d, 2*r, 'pca');
t2 = clock;
fprintf('get Ui cost %.1f seconds\n', etime(t2,t1));


n = size(X, 2);
display = 0;
if isfield(opts,'display'); display = opts.display; end

moveflag = true(1,n);
Mout = zeros(size(X));


for i = 1 : n
    x = X(:,i);
    t1 = clock;
    [Mout(:,i), moveflag(i)] = GD_PutM(x, data, Ui, d, r, 'xy', opts);
    t2 = clock;
    delta_t = etime(t2,t1);
    if display; fprintf('-----%d-th initial point costs %.2f seconds-----\n',i,delta_t); end
end

info.moveflag = moveflag;

end
    