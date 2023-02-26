import numpy as np
from noiseless_obj_grad import noiseless_obj_grad

def noiseless_PutM(P, Ui, x, d, r, algo, opts):
    
    # P : data collection
    # x : initial point to be moved to the manifold
    #
    # xiayq @ 8/19/2019
    #
    # xiayq@zju.edu.cn
    # refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228
    
    def getopts(opts):
        maxiter = 10; diff_tol = 0.01; display = False; alpha = 0.5; initer = 10
        if 'maxiter' in opts.keys():
            maxiter = opts['maxiter']
        if 'epsilon' in opts.keys():
            diff_tol = opts['diff_tol']
        if 'display' in opts.keys():
            display = opts['display']
        if 'alpha0' in opts.keys():
            alpha = opts['alpha0']
        if 'alpha' in opts.keys():
            alpha = opts['alpha']
        if 'initer' in opts.keys():
            initer = opts['initer']
        return [maxiter, diff_tol, display, alpha, initer]
    
    moveflag = True
    [maxiter, diff_tol, display, alpha0, initer] = getopts(opts)
    D = P.shape[0]
    
    monitor = False        
    for i in range(maxiter):
        
        if algo == 'GD':
            [flag, f_old, G, H] = noiseless_obj_grad(P, Ui, x, r)
            if not flag:
                moveflag = False
                monitor = True
        elif algo == 'SCGD':
            [flag, f_old, G, H] = noiseless_obj_grad(P, Ui, x, r)
            if not flag:
                moveflag = False
                monitor = True
            else:
                V = np.linalg.svd(H)[0]
                V = V[:, 0:(D-d)]
                G = V@(V.T@G)
        
        if monitor:
            pass
        else:
            if display:
                print('iter = %d: objective value is %.8f' % (i, f_old))
                    
            # gradient descent
            # G = G/np.linalg.norm(G)
            x_old = x.reshape(x.shape[0], 1); alpha = alpha0
            for iteration in range(10):
                x = x_old - alpha*G
                [flag, f, Matrix_G, Matrix_H] = noiseless_obj_grad(P, Ui, x, r)
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
                if f < 0:
                    raise ValueError('objective value is negative.')
                    
                if (abs(f - f_old) < diff_tol*max(f_old, np.finfo(np.float64).eps)):
                    break
    if monitor:
        pass
    else:
        if display:
            print('iter = %d: objective value is %.8f' % ((i + 1), f))
    
    return [x, moveflag]