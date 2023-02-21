function y = fun_ddh(x)

% second order derivative of h(x), where 
%          h(x) = exp(-1/t^2) for x>0 and 0 otherwise.
% the second order derivative should be
%        ddh(x) = (4t^{-6}-6t^{-4})exp(-1/t^2) for x>0 and 0 otherwise. 
%
% xiayq @ 8/14/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228




y = zeros(size(x));
idx = find(x > 0);
t = x(idx);
y(idx) = (4./t.^6 - 6./t.^4) .* exp(-1./t.^2);


% test code
% x1 = rand(1);
% delta = randn(1)*0.001;
% x2 = x1 + delta;
% 
% h1 = fun_dh(x1);
% h2 = fun_dh(x2);
% dh1 = fun_ddh(x1);
% 
% disp([h2, h1+dh1*delta])