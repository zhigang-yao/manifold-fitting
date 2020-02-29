function [P, d] = pdtorus(R, r, X)

% return the projection p of x onto the torus:
%    {(a,b,c) | (R-sqrt(a^2+b^2))^2+c^2 = r^2}
% d = ||p-x||_2
%
% xiayq @ 8/14/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228



temp1 = sqrt(X(1,:).^2+X(2,:).^2);
temp2 = R./temp1;
temp3 = (R-temp1).^2 + X(3,:).^2;
w1 = r./sqrt(temp3);
w0 = temp2+(1-temp2).*w1;

P = bsxfun(@times, X, [w0;w0;w1]);

d = sqrt(sum((X-P).^2));