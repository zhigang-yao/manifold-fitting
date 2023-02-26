import numpy as np
import time
from getUi import getUi
from noiseless_PutM import noiseless_PutM

def manfit_km17(data, d, r, X, opts):
    # This code is to implement the LPCA method in 
    # Mohammed, K. and Narayanan, H. (2017). Manifold Learning Using Kernel Density Estimation and Local Principal Components Analysis. 
    # arXiv:1709.03615.
    #
    # input : 
    # data - dataset
    # d    - dimension of manifold
    # r    - radius of neighborhood
    # X    - the set of initial points
    # opts - optional parameters, can be set as opts=[]
    #        maxiter : maxiter of GD
    #        epsilon : stopping criterian
    #        alpha0 : initial step length
    #        display : whether display iteration information or not
    #
    # output :
    # Mout - output manifold
    # info - informations of algorithm
    #
    # xiayq @ 8/18/2019 first version
    # xiayq @ 8/19/2019 change Pi to Ui, so that memory is in low level
    #
    # xiayq@zju.edu.cn
    # refered to Yao, Z. and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228

    t1 = time.time()
    Ui = getUi(data, data, d, r, 'pca')
    t2 = time.time()
    print('get Ui cost %.1f seconds' % (t2 - t1))
    
    n = X.shape[1]
    display = 0
    if 'display' in opts.keys():
        display = opts['display']
    
    moveflag = np.full_like(np.zeros((1, n)), True, dtype = bool)
    Mout = np.empty(((d + 1), 0))
    
    info = {}
    for i in range(n):
        x = X[:, i]
        t1 = time.time()
        [Mout_Part, moveflag[0, i]] = noiseless_PutM(data, Ui, x, d, r, 'SCGD', opts)
        Mout = np.hstack((Mout, Mout_Part))
        t2 = time.time()
        delta_t = t2 - t1
        if display:
            print('-----%d-th initial point costs %.2f seconds-----' % ((i + 1), delta_t))
    
    info['moveflag'] = moveflag
    
    return [Mout, info]