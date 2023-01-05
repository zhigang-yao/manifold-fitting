function y = fun_h(x)

% the function h(x)
%          h(x) = exp(-1/t^2) for x>0 and 0 otherwise.
%
% xiayq @ 8/14/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228



y = zeros(size(x));
idx = find(x > 0);
t = x(idx);
y(idx) = exp(-1./t.^2);

