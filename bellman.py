# -*- coding: utf-8 -*-

import numpy as np
from scipy import stats


def bin_segments(V_,max_k,criterion):
    V = stats.zscore(V_)
    T,I,J = np.unique(V,return_index=True,return_inverse=True)
    occurrences = np.histogram(J,max(J+1))
    n = len(T)
    cumocc = np.insert(np.cumsum(occurrences[0]),0,0)
    cs = np.insert(np.cumsum(np.multiply(T,occurrences[0])),0,0)
    css = np.insert(np.cumsum(np.multiply(np.power(T,2),occurrences[0])),0,0)
    sigma = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            sigma[i,j] = (css[j+1] - css[i]) - (1/(cumocc[j+1]-cumocc[i]))*np.power((cs[j+1]-cs[i]),2)   
    E = np.zeros((max_k,n))
    E[0,:] = sigma[0,:]
    B = np.zeros((max_k,n))
    cost = np.inf
    for k in range(1,max_k):
        for cutp in range(k,n):
            E[k,cutp] = np.min(np.add(E[k-1,0:cutp],sigma[1:cutp+1,cutp]))
            B[k,cutp] = int(np.argmin(np.add(E[k-1,0:cutp],sigma[1:cutp+1,cutp])))
        if criterion == "manual":
            c = E[k,n-1]
        elif criterion == "BIC":
            c = E[k,n-1] + k*np.log(len(V))
        elif criterion == "AIC":
            c = E[k,n-1] + 2*k
        if c < cost:
            opt_k = k
            cost = c
        else:
            opt_k = k-1
            break
    B = B.astype(int)
    cuts = [B[opt_k,n-1],n-1]
    for j in range(opt_k-1,0,-1):
        tmp = B[j,cuts[0]]
        cuts = np.insert(cuts,0,tmp)
    cuts = np.insert(np.unique(cuts),0,0)
    buckets = np.zeros(n)
    bounds = []
    U = np.unique(V_)
    for i in range(len(cuts)-1):
        buckets[cuts[i]+1:cuts[i+1]+1] = i
        bounds.append(U[cuts[i]+1]) 
    bounds = bounds[1:]
    assign = buckets[J]
        
    return assign, bounds, cost, opt_k
