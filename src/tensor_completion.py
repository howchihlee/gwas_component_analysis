import numpy as np
from scipy import linalg
def imode_cov(X, ind):
    # input: a tensor X and mode index ind
    # output: covariance matrix of X_(i) 
    Y = imode_fold(X, ind)
    return Y.dot(Y.T)  

def imode_unfold(X, ind):
    return np.swapaxes(X, 0, ind).reshape((X.shape[ind], -1))

def imode_fold(X, ind, dims):
    dims = list(dims)
    dims[ind], dims[0] = dims[0], dims[ind]
    return np.swapaxes(X.reshape(dims), ind, 0)

def imode_product(X, M, ind):
    # input X ndarray, M: matrix
    # return X\time_i M
    X1 = np.swapaxes(X, ind, X.ndim -1)
    M = np.reshape(M, [1] * (X.ndim-1) + [i for i in M.shape])
    Y = np.sum(np.expand_dims(X1, X.ndim) * M, axis = X1.ndim-1)
    return np.swapaxes(Y, ind, X.ndim-1)
        
def imode_soft_thresholding(X, lamb, ind):
    ## compute soft_thresholding on singular values of X
    Y = imode_unfold(X, ind)
    U, s, V = linalg.svd(Y, full_matrices=False)
    s = s - lamb
    s[s<0] = 0
    return imode_fold((s[np.newaxis, :] * U).dot(V), ind, X.shape)
    

def HaLRTC(T, ind_mat, alpha, rho, iters):
    X = T.copy()
    Y = np.zeros([X.ndim] + [i for i in X.shape])
    M = np.zeros([X.ndim] + [i for i in X.shape])
   
    for k in range(iters):
        rho *= 1.1
        for i in range(X.ndim):
            M[i, :] = imode_soft_thresholding(X + Y[i, :] / rho, alpha[i] / rho, i)
        X[ind_mat == 0] = np.mean(M - Y / rho, axis = 0)[ind_mat == 0]
        Y = Y - rho * (M - X[np.newaxis, :])
    return X
