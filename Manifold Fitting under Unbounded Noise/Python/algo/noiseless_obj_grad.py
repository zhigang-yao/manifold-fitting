import numpy as np
from fun_theta import fun_theta
from fun_dtheta import fun_dtheta
from fun_ddtheta import fun_ddtheta

def noiseless_obj_grad(P, Ui, x, r):
    # calculate the objective value f, that is, ||Pi_x F(x)||_2^2
    # and the gradient G of ||Pi_x F(x)||_2^2
    #
    # input variables
    # P : points for epsilon-net
    # Pi : Projection at each pi, a D*D*n tensor, Pi = I-(UiUi')
    # x : point to be projected to the manifold
    # d : dimension of intrinsic manifold
    # r : bandwidth parameter
    #
    # output variables
    # flag : whether x has enough neighbors or not
    # f : objective value
    # G : 0.5 times gradient
    #
    # variables in code
    # D : dimension of the ambient space
    # n : number of points in {pi}
    #
    # xiayq @ 8/19/2019 : change Pi to Ui to save memory
    #
    # xiayq@zju.edu.cn
    # refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228
    
    x = x.reshape((x.shape[0], 1))
    [D, n] = P.shape
    
    # calculate the difference bwtween x and each plane say bi=Pi^i(x-pi) 
    # and the squared distance fi = ||bi||_2^2, di=sqrt(fi)/2r
    bi = np.empty((D, 0))
    dxp = -(P - x)
    for i in range(n):
        bi_Part = (np.eye(D) - Ui[:, :, i]@Ui[:, :, i].T)@(x - P[:, i].reshape((P.shape[0], 1)))
        bi = np.hstack((bi, bi_Part))
    fi = sum(bi**2).reshape((1, bi.shape[1]))
    dists = sum(dxp**2).reshape((1, dxp.shape[1]))
    di = np.sqrt(dists)/(2*r)
    indexing = np.where(di <= 1)[1]
    num_idx = np.size(indexing)
    
    flag = (num_idx > 0)
    if flag == 0:
        f = -1
        G = np.zeros((D, 1))
        H = np.zeros((D, D))
        
    else:
        
        bi = bi[:, indexing]
        fi = fi[:, indexing]
        di = di[:, indexing]
        dists = dists[:, indexing]
        Ui = Ui[:, :, indexing]
        dxp = dxp[:, indexing]
        
        # calculate tilde_alpha_i, alpha_i, and alpha
        tilde_alpha_i = fun_theta(di)
        alpha = np.maximum(np.sum(tilde_alpha_i), np.finfo(np.float64).eps)
        alpha_i = tilde_alpha_i/alpha
        
        # calculate the objective value
        f = np.sum(fi*alpha_i)
        
        # calculate dalpha_i
        d_tilde_alpha_i = (1/(2*r)) / np.sqrt(dists) * fun_dtheta(di)
        d_tilde_alpha_i = dxp*d_tilde_alpha_i
        da = d_tilde_alpha_i.sum(axis = 1).reshape((d_tilde_alpha_i.shape[0], 1))
        dalpha_i = 1/alpha*(d_tilde_alpha_i - da*alpha_i)
        G = (dalpha_i*fi).sum(axis = 1) + 2*((bi*alpha_i).sum(axis = 1))
        G = G.reshape((G.shape[0], 1))
        
        # calculate Hessian of tilde_alpha_i
        # Hi = np.zeros((D, D, num_idx))
        Ax = np.zeros((D, D))
        sum_Hi = np.zeros((D, D))
        for i in range(num_idx):
            dp = dxp[:, i].reshape((dxp.shape[0], 1))
            a1 = -(dists[:, i]**(-1.5))*fun_dtheta(di[:, i])[0]
            a2 = 1/(2*r)/dists[:, i]*fun_ddtheta(di[:, i])[0]
            a3 = (dists[:, i]**(-0.5))*fun_dtheta(di[:, i])[0]
            # Hi[:, :, i] = ((a1 + a2)@(dp@dp.T) + a3)/(2*r)
            sum_Hi = sum_Hi + ((a1 + a2)*(dp@dp.T) + a3)/(2*r)
            Ax = Ax + alpha_i[:, i]*(np.eye(D) - Ui[:, :, i]@Ui[:, :, i].T)
        
        H = np.zeros((D, D))
        for i in range(num_idx):
            a = dalpha_i[:, i].reshape((dalpha_i.shape[0], 1)); b = bi[:, i].reshape((bi.shape[0], 1))
            dp = dxp[:, i].reshape((dxp.shape[0], 1))
            a1 = -(dists[:, i]**(-1.5))*fun_dtheta(di[:, i])[0]
            a2 = 1/(2*r)/dists[:, i]*fun_ddtheta(di[:, i])[0]
            a3 = (dists[:, i]**(-0.5))*fun_dtheta(di[:, i])[0]
            Hi = ((a1 + a2)*(dp@dp.T) + a3)/(2*r)
            H = H + 2*(a@b.T + b@a.T) - fi[:, i]/alpha*(da@a.T + a@da.T - Hi + alpha_i[:, i]*sum_Hi)
        H = H + 2*Ax
        # print(H)
    
    return [flag, f, G, H]