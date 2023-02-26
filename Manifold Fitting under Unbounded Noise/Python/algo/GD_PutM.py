from get_obj_grad import get_obj_grad
import time
import numpy as np

def GD_PutM(x, P, Ui, d, r, Ftype, opts):
    
    # each column of P is a point in X1, the epsilon-net
    # data are the sampled data, denoted by X0 in the paper
    #
    # xiayq @ 8/14/2019
    # xiayq@zju.edu.cn
    # refered to Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228
    
    def getopts(opts):
        maxiter = 10; epsilon = 0.01; display = False; alpha = 1; initer = 10
        if 'maxiter' in opts.keys():
            maxiter = opts['maxiter']
        if 'epsilon' in opts.keys():
            epsilon = opts['epsilon']
        if 'display' in opts.keys():
            display = opts['display']
        if 'alpha0' in opts.keys():
            alpha = opts['alpha0']
        if 'alpha' in opts.keys():
            alpha = opts['alpha']
        if 'initer' in opts.keys():
            initer = opts['initer']
        return [maxiter, epsilon, display, alpha, initer]
    
    moveflag = True
    [maxiter, epsilon, display, alpha0, initer] = getopts(opts)
    if 'beta' in opts.keys():
        beta = opts['beta']
    else:
        beta = d + 2
    
    monitor = False
    for i in range(maxiter):
        t1 = time.time()
        # [flag, f_old, G] = get_obj_grad_old(P, P, x, d, r, Ftype, 'pca')
        [flag, f_old, G] = get_obj_grad(x, P, Ui, d, r, Ftype, beta)
        t2 = time.time()
        delta_t = t2 - t1
        
        if not flag:
            moveflag = False
            x = x.reshape(x.shape[0], 1)
        else:
            if display:
                print('iter = %d: objective value is %.8f, costs %.2f seconds' % (i, f_old, delta_t))
            if 'logname' in opts.keys():
                with open(opts['logname'], 'a') as f:
                    f.write('iter = %d: objective value is %.8f, costs %.2f seconds' % (i, f_old, delta_t))
                    f.closed
                    
            # gradient descent
            # G = G/np.linalg.norm(G)
            x_old = x.reshape(x.shape[0], 1); alpha = alpha0
            t1 = time.time()
            for iteration in range(initer):
                x = x_old - alpha*G
                [flag, f, Matrix] = get_obj_grad(x, P, Ui, d, r, Ftype, beta)
                # [flag, f] = get_obj_grad_old(P, P, x, d, r, Ftype, 'pca')
                if ((flag == 0) | (f > f_old)) :
                    # if the movement of x is bad, decrease the step length or move
                    # back to x_old
                    if (iteration + 1) < initer:
                        alpha = alpha / 2
                    else:
                        x = x_old
                        monitor = True
                else:
                    break
            if monitor:
                pass
            else:
                t2 = time.time()
                delta_t = t2 - t1
                
                if f < 0:
                    raise ValueError('objective value is negative.')
                    
                if ((i > 1) & (abs(f) < np.finfo(np.float64).eps)):
                    break
                # if (abs(f-f_old) < np.finfo(np.float64).eps*f_old):
                #     break
    if monitor:
        pass
    else:
        if display:
            print('iter = %d: objective value is %.8f, costs %.2f seconds' % ((i + 1), f, delta_t))
        if 'logname' in opts.keys():
            with open(opts['logname'], 'a') as f:
                f.write('iter = %d: objective value is %.8f, costs %.2f seconds' % ((i + 1), f, delta_t))
                f.closed
    
    return [x, moveflag]