import numpy as np
import scipy.linalg

def spec_thresholding_fat(X, lamb):
    cov = X.dot(X.T)
    #cov = reduce(lambda x, y: x+y, [np.zeros_like(covs[0])] + covs)

    D, U = scipy.linalg.eigh(cov)
    s = np.sqrt(np.abs(D))
    Dhat = (s - lamb) / (s + 1e-6)
    Dhat[s < lamb] = 0
    T = (Dhat[np.newaxis, :] * U).dot(U.T)
    Xt = T.dot(X)
    
    return Xt

def spec_thresholding(X, lamb):
    U, s, V = scipy.linalg.svd(X, full_matrices=False)
    s = np.maximum(s - lamb, 0)
    return np.dot(U, np.dot(np.diag(s), V))

def matrix_imputing(Y, ind, lamb, max_iter):
    B = np.zeros(Y.shape)
    iter = 0
    while iter < max_iter:
        U, s, V = np.linalg.svd(ind * Y + (1-ind) * B, full_matrices=False)
        s = s - lamb
        s[s<0] = 0
        B = np.dot(U, np.dot(np.diag(s), V))
        iter += 1
    return B

def matrix_imputing_ADMM(ind, M, beta, X0, epsilon, max_iter, compute_obj = False, fat = False):
    ## algorhtm 2 in Hu et al. 2013, Fast and accurate matrix completion 
    ## via truncated nuclear norm regulatization, 
    Xk = X0.copy()
    Wk = Xk.copy()
    Yk = Xk.copy()
    Xold = np.zeros(Xk.shape)
    iters = 0
    obj_fun = []
    t = 1.1   
    
    while iters < max_iter:
        Xold = Xk.copy()
        if fat:
            Xk = spec_thresholding_fat(Wk - Yk / beta, 1.0 / beta)
        else:
            Xk = spec_thresholding(Wk - Yk / beta, 1.0 / beta)
        
        Wk = Xk +  Yk / beta
        Wk[ind == 1] = M[ind == 1]
        Yk += (Xk - Wk) * beta
        
        if compute_obj:
            U, s, V = scipy.linalg.svd(Xk, full_matrices=False)
            tmp = np.sum(s)
            #tmp += np.linalg.norm(Xk - Wk) * beta / 2.0 + np.sum(Yk * (Xk - Wk))
            obj_fun.append(tmp)
        iters += 1
        beta *= t
    return Xk, obj_fun
