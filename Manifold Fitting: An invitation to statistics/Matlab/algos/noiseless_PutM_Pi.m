function [x, moveflag] = noiseless_PutM_Pi(P, Pi, x, d, r, algo, opts)

% P : data collection
% x : initial point to be moved to the manifold
%
% xiayq @ 8/19/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228


moveflag = true;
[maxiter, diff_tol, display, alpha0, initer] = getopts(opts);
D = size(P,1);

for i  =  1 : maxiter
    
    switch algo
        case 'GD'
            [flag, f_old, G] = noiseless_obj_grad_Pi(P, Pi, x, r);
            if ~flag; moveflag = false; return; end
        case 'SCGD'
            [flag, f_old, G, H] = noiseless_obj_grad_Pi(P, Pi, x, r);
            if ~flag; moveflag = false; return; end
            [V, ~, ~] = svd(H);
            V = V(:, 1:D-d);
            G = V*(V'*G);
    end
    if display; fprintf('iter = %d: objective value is %.8f\n', i-1, f_old); end
   
    % gradient descent
    %G = G/norm(G,2); 
    x_old = x; alpha = alpha0;
    for iter = 1 : 10
        x = x_old - alpha*G;
        [flag, f] = noiseless_obj_grad_Pi(P, Pi, x, r);
        if flag == 0 || f > f_old 
            % if the movement of x is bad, decrease the step length or move
            % back to x_old
            if iter < initer
                alpha = alpha / 2;
            else
                x = x_old;
                return;
            end
        else
            break;
        end
    end
    
    if f < 0
        error('objective value is negative.');
    end
   
    if abs(f-f_old) < diff_tol*max(f_old,eps); break; end
    
end

 if display; fprintf('iter = %d: objective value is %.8f\n', i, f); end

end

function [maxiter, diff_tol, display, alpha, initer] = getopts(opts)
    maxiter = 10; diff_tol = 1e-2; display = false; alpha = 0.5; initer = 10;
    if isfield(opts, 'maxiter'); maxiter = opts.maxiter; end
    if isfield(opts, 'epsilon'); diff_tol = opts.diff_tol; end
    if isfield(opts, 'display'); display = opts.display; end
    if isfield(opts, 'alpha0'); alpha = opts.alpha0; end
    if isfield(opts, 'alpha'); alpha = opts.alpha; end
    if isfield(opts, 'initer'); initer = opts.initer; end
    
end