import numpy as np

def get_obj_grad(x, P, Ui, d, r, Ftype, beta):
    # calculate the objective value f, that is, ||Pi_x F(x)||_2^2
    # and the gradient G of ||Pi_x F(x)||_2^2
    #
    # input variables
    # P : points for epsilon-net
    # data : noisy samples
    # x : point to be projected to the manifold
    # d : dimension of intrinsic manifold
    # r : bandwidth parameter
    # Ftype : Fi is selected as Fi = x-pi ('xy') or Fi = Pi(x-pi) ('cf')
    #
    # output variables
    # flag : whether x has enough neighbors or not
    # f : objective value
    # G : 0.5 times gradient
    #
    # variables in code
    # D : dimension of the ambient space
    #
    # xiayq @ 8/14/2019
    # xiayq@zju.edu.cn
    # refered to Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228
    
    # def Ax_Ui(Ui, alpha_i, v):
    #     Ui = np.array(Ui)
    #     [D, n] = Ui.shape[[0, 2]]
    #     Axv = np.zeros((D, 1))
    #     for i in range(n):
    #         ui = Ui[:, :, i]
    #         Axv = Axv + (v - ui@(ui.T@v))*alpha_i[i]
    #     return Axv
    
    x = x.reshape(x.shape[0], 1)
    D = P.shape[0]
    
    # calculate x - pi
    dxp = -(P - x)
    dist2 = sum(dxp**2).reshape(1, dxp.shape[1])
    
    indexing = np.where(dist2 < (r*r))[1]
    num_idx = np.size(indexing)
    
    flag = (num_idx > 0)
    if flag == 0:
        f = -1; G = np.zeros((D, 1))
    
    else:
        # P = P[:, indexing]
        dxp = dxp[:, indexing]
        dist2 = dist2[:, indexing].reshape(num_idx, )
        Ui = Ui[:, :, indexing]
        
        # calculate tilde_alpha_i, alpha_i and alpha
        tilde_alpha_i = np.zeros((num_idx, 1))
        for i in range(num_idx):
            tilde_alpha_i[i] = (1 - dist2[i]/(r*r))**beta
        alpha = sum(tilde_alpha_i)
        alpha_i = tilde_alpha_i/alpha
        
        # calculate Pi
        # Pi = np.zeros((D, D, num_idx))
        # for i in range(num_idx):
        #     Pi[:, :, i] = getPi(P[:, i], data, d, r, dr_type)
        
        # calculate Ax and Fi = alpha_i Pi(x-pi)
        Fi = np.zeros((D, num_idx))
        Ax = np.zeros(D)
        for i in range(num_idx):
            Ax = Ax + alpha_i[i]*(np.eye(D) - Ui[:, :, i]@Ui[:, :, i].T)
            if Ftype == 'cf':
                Pi = np.eye(D)- Ui[:, :, i]@Ui[:, :, i].T
                Fi[:, i] = Pi@dxp[:, i]*alpha_i[i]
                # Fi[:, i] = (dxp[:, i] - Ui[:, :, i]@(Ui[:, :, i].T@dxp[:, i]))*alpha_i[i]
            elif Ftype == 'xy':
                Fi[:, i] = dxp[:, i]*alpha_i[i]
        Ax = (Ax + Ax.T)/2
        
        # calculate Fx
        Fx = Fi.sum(axis = 1).reshape(Fi.shape[0], 1)
        
        # calculate V and Pix*b
        l, V = np.linalg.eig(Ax)
        [l, idx] = [l[np.argsort(-l)], np.argsort(-l)]
        V = V[:, idx]
        
        # Ps = np.zeros((D, D, D))
        # for i in range(D)
        #     Ps[:, :, i] = V[:, i]@V[:, i].T
        # Pix = Ps[:, :, 0:D-d].sum(axis = 2)
        
        # calulate the objective value
        # b = Pix@Fx
        b = (Fx - V[:, (D-d):]@(V[:, (D-d):].T@Fx)).reshape(Fx.shape[0], 1)
        f = sum(b**2)
        
        # calculate B and T
        B = b@Fx.T
        
        T = np.zeros(D)
        v = V[:, 0:(D-d)]
        for j in range((D-d), D):
            vj = V[:, j].reshape(V.shape[0], 1)
            dlambda = 1/(l[j] - l[0:D-d])
            temp = v*dlambda.T@v.T
            T = T + temp@(B@vj)@vj.T + vj@(vj.T@B)@temp
            
        # for i in range(D-d):
        #     for j in range(D-d, D):
        #         vi = V[:, i]
        #         vj = V[:, j]
        #         T = T + 1/(l[j] - l[i])*(vi@(vi.T@B@vj)@vj.T + vj@(vj.T@B@vi)@vi.T))
        
        # calculate dalpha_i
        bar_alphai = tilde_alpha_i**((beta-1)/beta)
        dalpha_i = dxp*bar_alphai.T
        dalpha_i = dalpha_i - dalpha_i.sum(axis = 1).reshape(dalpha_i.shape[0], 1)*alpha_i.T
        dalpha_i = -dalpha_i*(2*beta/(r*r*alpha))
        
        # calculate <T, Pi>, <b, Fi>
        PiT = np.zeros((1, num_idx))
        for i in range(num_idx):
            # Pi = np.eye(D) - Ui[:, :, i]@Ui[:, :, i].T
            # PiT[i] = sum(Pi[:]*T[:])
            ui = Ui[:, :, i]
            Tui = T@ui
            PiT[:, i] = T.trace() - sum(ui.flatten()*Tui.flatten())
        Fib = b.T@Fi
        
        if Ftype == 'cf':
            G = (dalpha_i*(Fib + PiT)).sum(axis = 1).reshape(dalpha_i.shape[0], 1) + Ax@b
        elif Ftype == 'xy':
            G = (dalpha_i*(Fib + PiT)).sum(axis = 1).reshape(dalpha_i.shape[0], 1) + b
        
    return [flag, f, G]