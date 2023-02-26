import math
import numpy as np
import pickle
import time
from manfit_our import manfit_our
from manfit_cf18 import manfit_cf18
from manfit_km17 import manfit_km17
import matplotlib.pyplot as plt
import scipy.io

# xiayq @ 12/11/2022
# xiayq0121@zufe.edu.cn
# refered to Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228

# parameters for data
D = 3; Num = 1000
NumTrials = 20
tau = 1; sigma = 0.02 #0.04

# method setup
algos = ['ours', 'cf18', 'km17']

if sigma == 0.09:
    r1 = 2*math.sqrt(sigma)
    r2 = 2.5*math.sqrt(sigma) 
    r3 = 1.2*math.sqrt(sigma)
elif sigma == 0.08:
    r1 = 2*math.sqrt(sigma) 
    r2 = 2.8*math.sqrt(sigma)
    r3 = 1.2*math.sqrt(sigma)
elif sigma == 0.07:
    r1 = 2*math.sqrt(sigma)
    r2 = 2*math.sqrt(sigma)
    r3 = 1.2*math.sqrt(sigma)
elif sigma == 0.06:
    r1 = 2*math.sqrt(sigma)
    r2 = 2.5*math.sqrt(sigma)
    r3 = 1.2*math.sqrt(sigma)
elif sigma == 0.05:
    r1 = 2*math.sqrt(sigma)
    r2 = 2.5*math.sqrt(sigma)
    r3 = 1.2*math.sqrt(sigma)
elif sigma == 0.04:
    r1 = 2*math.sqrt(sigma)
    r2 = 3*math.sqrt(sigma)
    r3 = 1.2*math.sqrt(sigma)
elif sigma == 0.03:
    r1 = 2.2*math.sqrt(sigma)
    r2 = 2.8*math.sqrt(sigma)
    r3 = 1.4*math.sqrt(sigma)
elif sigma == 0.02:
    r1 = 2.5*math.sqrt(sigma)
    r2 = 2.5*math.sqrt(sigma)
    r3 = 1.6*math.sqrt(sigma)
elif sigma == 0.01:
    r1 = 3*math.sqrt(sigma)
    r2 = 3.6*math.sqrt(sigma)
    r3 = 2*math.sqrt(sigma)
    
num_algo = np.size(algos)

avgdists = -np.ones((num_algo, NumTrials))
maxdists = -np.ones((num_algo, NumTrials))
ts = -np.ones((num_algo, NumTrials))

Mouts = np.empty((num_algo,NumTrials), dtype = object)
proj_Mouts = np.empty((num_algo,NumTrials), dtype = object)
infos = np.empty((num_algo,NumTrials), dtype = object)
Dist2 = np.empty((num_algo,NumTrials), dtype = object)
Dist2_move = np.empty((num_algo,NumTrials), dtype = object)

for rep in range(NumTrials):
    
    print('------ Trial %d ------' % (rep + 1))
    
    fname = 'simulations\\sphere\\t%d_n%d_s%.2f_trial%d.pkl' % (tau, Num, sigma, (rep + 1))
    
    # Loading data from MATLAB
    # Data = scipy.io.loadmat(fname)
    # samples = Data['samples']
    # data_ini = Data['data_ini']
    
    try:
        with open (fname, 'rb') as f:
            Data = pickle.load(f)
        NumSample = Data['samples'].shape[1]
        NumIni = Data['data_ini'].shape[1]
        print('load data %d' % (rep + 1))
        
    except:
        
        # generate data
        NumSample = Num
        NumIni = Num
        np.random.seed(rep + 1)
        samples = np.random.randn(3, NumSample)
        samples = samples@np.diag(1/np.sqrt(sum(samples**2)))+ sigma*np.random.randn(3, NumSample)
        
        data_ini = np.random.randn(3, NumIni)
        data_ini = data_ini@np.diag(1/np.sqrt(sum(data_ini**2)))+0.5*math.sqrt(sigma)/math.sqrt(D)*(2*np.random.rand(3, NumIni) - 1)
        
        with open (fname, 'wb') as f:
            pickle.dump({'samples': samples, 'data_ini': data_ini}, f)
        print('generate data %d' % (rep + 1))
    
    # parameters for algorithm
    dim = 2
    opts = {}
    opts['epsilon'] = 1e-16
    opts['diff_tol'] = 1e-4
    opts['maxiter'] = 50
    opts['display'] = 0
    opts['initer'] = 10
    
    for i in range(num_algo):
        algo = algos[i]
        
        t1 = time.time()
        if algo == 'ours':
            [Mout, info] = manfit_our(samples, dim, r1, data_ini, opts)
        elif algo == 'cf18':
            [Mout, info] = manfit_cf18(samples, dim, r2, data_ini, opts)
        elif algo == 'km17':
            [Mout, info] = manfit_km17(samples, dim, r3, data_ini, opts)
        # elif algo == 'uo11':
            # Mout = pc_project_multidim(samples, data_ini, r, dim)
            # Mout = Mout.T
            # info['moveflag'] = np.full_like(np.zeros((1, Mout.shape[1])), True)
        t2 = time.time()
        
        Mouts[i, rep] = Mout
        infos[i, rep] = info
        ts[i, rep] = t2 - t1
        
        print('Trial %d with algo %s costs %.2f seconds' % ((rep + 1), algo, (t2 - t1)))
        
    # fname = 'out/sphere/Dist_t%d_s%.2f_r%.2f.pkl' % (tau, sigma, r)
    # with open (fname, 'wb') as f:
    #     pickle.dump({'Mouts': Mouts, 'infos': infos, 'ts': ts, 'Dist2': Dist2,
    #                  'Dist2_move': Dist2_move, 'avgdists': avgdists,
    #                  'maxdists': maxdists}, f)
    
    # calculate distances
    for i in range(num_algo):
        Mout = Mouts[i, rep]
        moveflag = infos[i, rep]['moveflag']
        
        proj_Mout = Mout*tau/np.sqrt(sum(Mout**2))
        proj_Mouts[i, rep] = proj_Mout
        Dist2[i, rep] = np.sqrt(sum((Mout - proj_Mout)**2))
        temp = Dist2[i, rep][moveflag.reshape((moveflag.shape[1],))]
        Dist2_move[i, rep] = temp
        
        print('Trial %d, algo %s: Max = %.8f, Avg = %.8f, Num = %d' % ((rep + 1), algos[i], np.max(temp), np.mean(temp), np.size(temp)))
        
        avgdists[i, rep] = np.mean(temp)
        maxdists[i, rep] = np.max(temp)
        
fname = 'out\\sphere\\Dist_n%d_s%.2f.pkl' % (Num, sigma)
with open (fname, 'wb') as f:
    pickle.dump({'algos': algos, 'Mouts': Mouts, 'proj_Mouts': proj_Mouts, 'infos': infos,
                 'ts': ts, 'Dist2': Dist2, 'Dist2_move': Dist2_move, 'avgdists': avgdists,
                 'maxdists': maxdists, 'r1': r1, 'r2': r2, 'r3': r3}, f)
    
# boxplot
MatName = 'out\\sphere\\Dist_n%d_s%.2f.pkl' % (Num, sigma)
PyName = 'out\\sphere\\Manapprox_Dist_n%d_s%.2f.mat' % (Num, sigma)

with open (MatName, 'rb') as f:
    D1 = pickle.load(f)
D2 = scipy.io.loadmat(PyName)
    
algos = D1['algos']
algos.append('ya21(deg=1)')
algos.append('ya21(deg=2)')
num_algo = np.size(algos)

maxdists = np.vstack((D1['maxdists'], D2['max_dist']))
NumTrials = maxdists.shape[1]

ordered = [4, 0, 3, 1, 2]
    
if NumTrials > 1:
    
    fig, ax = plt.subplots()
    plt.boxplot(maxdists[ordered,:].T)
    ax.set_xticklabels(np.array(algos)[ordered])
    plt.title('d = ' + str(dim) + ', ' + chr(963) + ' = ' + str(sigma) + ', n = ' + str(Num))
    fig.savefig('figures\\sphere_max_n%d_s%.2f.png' % (Num, sigma))
    
# Check the result differences between Python and MATLAB
# New = 'out\\sphere\\Dist_n%d_s%.2f.pkl' % (Num, sigma)
# Origin = 'out\\sphere\\Dist_n%d_s%.2f.mat' % (Num, sigma)

# with open (New, 'rb') as f:
#     New_Results = pickle.load(f)
# Origin_Results = scipy.io.loadmat(Origin)

# algos = ['ours', 'cf18', 'km17']
# num_algo = np.size(algos)

# Differences = np.empty((num_algo, NumTrials), dtype = object)
# for i in range(num_algo):
#     for rep in range(NumTrials):
#         Difference = New_Results['Dist2_move'][i, rep] - Origin_Results['Dist2_move'][i, rep]
#         Differences[i, rep] = Difference
        
# fname = 'out\\sphere\\Differences_n%d_s%.2f.pkl' % (Num, sigma)
# with open (fname, 'wb') as f:
#     pickle.dump(Differences, f)