function y = fun_dtheta(x)

% the first order deriative of the bump function
%       theta(x) = h(1-x)./(h(1-x)+h(x-1/4))
% that is
%       theta(x) = -h'(1-x)./(h(1-x)+h(x-1/4))
%                - h(1-x)(h'(x-1/4)-h'(1-x))./(h(1-x)+h(x-1/4)).^2
%
% xiayq @ 8/14/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228



x1 = 1-x; x2 = x-1/4;
h1 = fun_h(x1); h2 = fun_h(x2); h = h1 + h2;
dh1 = fun_dh(x1); dh2 = fun_dh(x2);

y = -dh1./h - h1.*(dh2-dh1)./(h.^2);

% test code
% x1 = rand(1)*0.3+0.4;
% delta = randn(1)*0.001;
% x2 = x1 + delta;
% 
% h1 = fun_theta(x1);
% h2 = fun_theta(x2);
% dh1 = fun_dtheta(x1);
% 
% disp([h2, h1+dh1*delta])


