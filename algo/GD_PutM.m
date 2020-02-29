function [x, moveflag] = GD_PutM(x, P, Ui, d, r, Ftype, opts)

% each column of P is a point in X1, the epsilon-net
% data are the sampled data, denoted by X0 in the paper
%
% xiayq @ 8/14/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228


moveflag = true;
[maxiter, epsilon, display, alpha0, initer] = getopts(opts);
if isfield(opts,'beta')
    beta = opts.beta;
else
    beta = d+2;
end

if isfield(opts,'logname'); fp = fopen(opts.logname,'a');end

for i  =  1 : maxiter
    t1 = clock;
    %[flag, f_old, G] = get_obj_grad_old(P, P, x, d, r, Ftype, 'pca');
    [flag, f_old, G] = get_obj_grad(x, P, Ui, d, r, Ftype, beta);
    t2 = clock;
    delta_t = etime(t2,t1);

    
    if ~flag; moveflag = false; return; end
    if display; fprintf('iter = %d: objective value is %.8f, costs %.2f seconds \n', i-1, f_old, delta_t); end
    if isfield(opts,'logname'); fprintf(fp,'iter = %d: objective value is %.8f, costs %.2f seconds \n', i-1, f_old, delta_t); end
        
   
    % gradient descent
    %G = G/norm(G,2); 
    x_old = x; alpha = alpha0;
    t1 = clock;
    for iter = 1 : initer
        x = x_old - alpha*G;
        [flag, f] = get_obj_grad(x, P, Ui, d, r, Ftype, beta);
        %[flag, f] = get_obj_grad_old(P, P, x, d, r, Ftype, 'pca');
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
    t2 = clock;
    delta_t = etime(t2,t1);
    
    if f < 0
        error('objective value is negative.');
    end
   
    if i > 1 && abs(f) < epsilon; break; end
    %if abs(f-f_old) < epsilon*f_old; break; end
    
end

 if display; fprintf('iter = %d: objective value is %.8f, costs %.2f seconds\n', i, f, delta_t); end
 if isfield(opts,'logname')
     fprintf(fp,'iter = %d: objective value is %.8f, costs %.2f seconds\n', i, f, delta_t); 
     fclose(fp);
 end
    

end

function [maxiter, epsilon, display, alpha, initer] = getopts(opts)
    maxiter = 10; epsilon = 1e-2; display = false; alpha = 1; initer=10;
    if isfield(opts, 'maxiter'); maxiter = opts.maxiter; end
    if isfield(opts, 'epsilon'); epsilon = opts.epsilon; end
    if isfield(opts, 'display'); display = opts.display; end
    if isfield(opts, 'alpha0'); alpha = opts.alpha0; end
    if isfield(opts, 'alpha'); alpha = opts.alpha; end
    if isfield(opts, 'initer'); initer = opts.initer; end
end