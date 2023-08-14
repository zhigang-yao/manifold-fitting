function y = fun_ddtheta(x)

% the second order deriative of the bump function
%       theta(x) = h(1-x)./(h(1-x)+h(x-1/4))
% that is
%       theta(x) = h''(1-x)./(h(1-x)+h(x-1/4)) 
%                - 2h'(1-t)(h'(1-t)-h'(t-1/4))./ (h(1-x)+h(x-1/4)).^2
%                - h(1-x)(h''(x-1/4)+h''(1-x))./(h(1-x)+h(x-1/4)).^2
%                + 2h(1-x)(h'(t-1/4)-h'(1-t)).^2./(h(1-x)+h(x-1/4)).^3
%
% xiayq @ 8/14/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228



x1 = 1-x; x2 = x-1/4;
h1 = fun_h(x1); h2 = fun_h(x2); h = h1 + h2;
dh1 = fun_dh(x1); dh2 = fun_dh(x2); dh = dh1-dh2;
ddh1 = fun_ddh(x1); ddh2 = fun_ddh(x2);

y = ddh1./h - (2*dh1.*dh + h1.*(ddh2+ddh1))./h.^2 + 2.*h1.*(dh.^2)./(h.^3);

% test code
% x1 = rand(1)*0.3+0.4;
% delta = randn(1)*0.0001;
% x2 = x1 + delta;
% 
% h1 = fun_dtheta(x1);
% h2 = fun_dtheta(x2);
% dh1 = fun_ddtheta(x1);
% 
% disp([h2, h1+dh1*delta])