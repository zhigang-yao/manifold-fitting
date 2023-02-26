import numpy as np

def pdtorus(R, r, X):
    
    # return the projection p of x onto the torus:
    #    {(a,b,c) | (R-sqrt(a^2+b^2))^2+c^2 = r^2}
    # d = ||p-x||_2
    #
    # xiayq @ 8/14/2019
    #
    # xiayq@zju.edu.cn
    # refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228
    
    temp1 = np.sqrt(X[0, :]**2 + X[1, :]**2)
    temp2 = R/temp1
    temp3 = (R - temp1)**2 + X[2, :]**2
    w1 = r/np.sqrt(temp3)
    w0 = temp2 + (1 - temp2)*w1
    
    P = X*np.vstack((w0, w0, w1))
    
    d = np.sqrt(sum((X - P)**2))
    
    return [P, d]