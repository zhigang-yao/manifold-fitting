import numpy as np
import math

def torusUnif(n, R, r):
    # generate n samples uniformly from a torus
    # R : the radius from the center of the hole to the center of the torus tube
    # r : the radius of the torus tube. 
    #
    # There is an R package with the same name
    # Algorithm 1 of Diaconis P, Holmes S, and Shahshahani M (2013). "Sampling from a manifold." 
    # Advances in Modern Statistical Theory and Applications: A Festschrift in honor of Morris L. Eaton. 
    # Institute of Mathematical Statistics, 102-125.
    
    count = 0
    theta = -np.ones((1, n))
    
    while count < n:
        xvec = np.random.rand(1)*2*math.pi
        yvec = np.random.rand(1)/math.pi
        fx = (1 + r/R*math.cos(xvec))/(2*math.pi)
        if yvec < fx:
            theta[:, count] = xvec
            count = count + 1
    
    phi = np.random.rand(1, n)*2*math.pi
    x = (R + r*np.cos(theta))*np.cos(phi)
    y = (R + r*np.cos(theta))*np.sin(phi)
    z = r*np.sin(theta)

    X = np.vstack((x, y, z))
    
    return X