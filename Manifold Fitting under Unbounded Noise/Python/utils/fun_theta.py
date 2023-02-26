import numpy as np
from fun_h import fun_h
# import matplotlib.pyplot as plt

def fun_theta(x):
    
    # the bump function theta satisfying 
    #       theta(x)=1 for x<1/4 and theta(x) = 0 for x \geq 1
    #       theta(x) = h(1-x)./(h(1-x)+h(x-1/4))
    #
    # xiayq @ 8/14/2019
    #
    # xiayq@zju.edu.cn
    # refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228
    
    x = np.array(x)
    t = fun_h(1 - x)
    y = t / (t + fun_h(x - 1/4))
    
    return y

# test code
# x = np.arange(0, 2, 0.01)
# y = fun_theta(x)
# plt.plot(x, y)