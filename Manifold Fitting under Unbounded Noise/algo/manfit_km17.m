function [Mout, info] = manfit_km17(data, d, r, X, opts)
% This code is to implement the LPCA method in 
% Mohammed, K. and Narayanan, H. (2017). Manifold Learning Using Kernel Density Estimation and Local Principal Components Analysis. 
% arXiv:1709.03615.
%
% input : 
% data - dataset
% d    - dimension of manifold
% r    - radius of neighborhood
% X    - the set of initial points
% opts - optional parameters, can be set as opts=[]
%        maxiter : maxiter of GD
%        epsilon : stopping criterian
%        alpha0 : initial step length
%        display : whether display iteration information or not
%
% output :
% Mout - output manifold
% info - informations of algorithm
%
% xiayq @ 8/18/2019 first version
% xiayq @ 8/19/2019 change Pi to Ui, so that memory is in low level
%
% xiayq@zju.edu.cn
% refered to Yao, Z. and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228

t1 = clock;
Ui = getUi(data, data, d, r, 'pca');
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
    [Mout(:,i), moveflag(i)] = noiseless_PutM(data, Ui, x, d, r, 'SCGD', opts);
    t2 = clock;
    delta_t = etime(t2,t1);
    if display; fprintf('-----%d-th initial point costs %.2f seconds-----\n',i,delta_t); end
end

info.moveflag = moveflag;

end
    