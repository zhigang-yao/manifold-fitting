def manfit_ours(sample, sig, sample_init, op_average=1):
    import numpy as np
    from scipy.spatial.distance import pdist, squareform
    from sklearn.neighbors import NearestNeighbors

    Mout = np.copy(sample_init)
    N = sample_init.shape[0]
    N0 = sample.shape[0]
    ns = np.arange(N0)

    r = 5 * sig / np.log10(N0)
    R = 10 * sig * np.sqrt(np.log(1 / sig)) / np.log10(N0)

    for ii in range(N):
      
        x = sample_init[ii, :]

        dists = squareform(pdist(np.vstack([x, sample])))[0, 1:]

        IDX1 = dists < 2 * r
        IDX1 = ns[IDX1]

        nbrs = NearestNeighbors(n_neighbors=5).fit(sample)
        IDX2 = nbrs.kneighbors(x.reshape(1, -1), return_distance=False).flatten()

        IDX = np.union1d(IDX1, IDX2)

        BNbr = sample[IDX, :]

        xbar = np.mean(BNbr, axis=0) + np.finfo(float).eps

        dx = x - xbar
        dx = dx / np.linalg.norm(dx)

        Q = np.linalg.qr(np.column_stack([dx, np.eye(dx.size)]))[0]

        sample_s = sample - x
        sample_s = sample_s @ Q

        CNbr = (np.abs(sample_s[:, 0]) < R) & (np.sum(sample_s[:, 1:] ** 2, axis=1) < r ** 2)

        if np.sum(CNbr) > 10:
            Mout[ii, :] = np.mean(sample[CNbr, :], axis=0)
        else:
            Mout[ii, :] = xbar

    return Mout

