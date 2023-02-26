import numpy as np
from scipy import linalg

def getUi(P, data, d, r, dr_type):
    # get Pi = UiUi' for each pi
    # each column of P is a point in X1, the epsilon-net(cf18)
    # data is the entire dataset
    #
    # xiayq @ 04/21/2022
    #       modify the comments
    #       remove 'cf1', replace 'cf2' by 'cf'
    # xiayq @ 08/18/2019 
    #       first version
    # xiayq0121@zufe.edu.cn
    # refered to Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228
    
    def selectNeighbor(data, p, d):
        N = data.shape[1]
        if N <= d:
            P0 = data
        else:
            idx = np.zeros((d, ), dtype = int)
            delta = abs(1 - np.sqrt(sum((data - p)**2)))
            M = np.zeros((d, N))
            M[0, :] = delta
            for i in range(d):
                curID = np.argmin(np.max(M, axis = 0))
                idx[i] = curID
                M[i, curID] = float('inf')
                x = data[:, curID]; x = x/np.linalg.norm(x)
                if i < (d - 1):
                    M[(i + 1), :] = abs(sum(x.reshape((x.shape[0], 1))*data))
            P0 = data[:, idx]
        return P0
    
    [n, N] = P.shape
    Ui = np.zeros((n, d, N))
    
    # rescale data with 1/r
    P = P/r; data = data/r
    
    for i in range(N):
        p = P[:, i].reshape((n, 1))
        dis2 = sum((data - p)**2)
        [dis2, idx] = [np.sort(dis2), np.argsort(dis2)]
        
        index = next(index for index, dis in enumerate(dis2) if dis > 1)
        
        if dis2[0] < np.finfo(np.float64).eps: # delete the point p itself
            # print('delete p%d from its neighbor' % (i + 1))
            P0 = data[:, idx[1:max(index, d+1)]]
        else:
            P0 = data[:, idx[0:max(index, d)]]
            
        if dr_type == 'pca':
            U = linalg.svd((P0 - p))[0]
            U = U[:, 0:d]
        elif dr_type == 'cf':
            # P0 - p takes p as the origin
            P0 = selectNeighbor((P0 - p), np.zeros(p.shape), d)
            U = linalg.orth(P0)
        # elif dr_type == 'cf1':
            # P0 = selectNeighbor(P0, p, d)
            # U = linalg.orth((P0 - p))
        Ui[:, :, i] = U
        
    return Ui