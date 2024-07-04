import numpy as np

def covis_covar_hamming(a,b,r,window,overlap):
    [N,M]=a.shape
    ham_win = np.array([np.hamming(window)]).T@np.ones((1,M))/np.sum(np.hamming(window))
    nbins = int(np.floor(N/(window-overlap))-1)    # Last range bin is removed
    cab = a*np.conj(b)
    cov = np.zeros((nbins,M))
    E1 = np.zeros((nbins,M))
    E2 = np.zeros((nbins,M))
    rc = np.zeros((nbins,1))

    for m in range(nbins):
        j = int(m*(window-overlap) + 1)
        k = int(j + window - 1)
        if k > N:
            k = N
        cabjk = cab[j-1:k,:]
        cov[m,:] = np.sum(cabjk*ham_win,axis=0)
        E1[m,:] = np.sum(np.power(np.abs(a[j-1:k,:]),2)*ham_win,axis=0)
        E2[m,:] = np.sum(np.power(np.abs(b[j-1:k,:]),2)*ham_win,axis=0)
        rc[m] = np.mean(r[j-1:k])
    return cov, E1, E2, rc