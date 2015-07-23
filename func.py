'''
Created on 02.07.2012

@author: Marina
'''
import numpy as np
import math
import utils


def extend_callfunc(x, Q, P, M, L, r, mu, sigma, sm, maxdegree, vecsize, experiment, opt_var, small_k):
    if (opt_var[0] == 1):
        r = utils.converter(x[:vecsize[0]], P, M)
    if (opt_var[1] == 1):
        sm = utils.listconverter(x[vecsize[0] :vecsize[1]], P, M)
    if (opt_var[2] == 1):
        sigma = utils.listconverter(x[vecsize[1]:vecsize[2]], P, M)
    if (opt_var[3] == 1):
        mu = utils.listconverter(x[vecsize[2]:], P, M)
    fval = extend_func_POIM(Q, P, M, L, r, mu, sigma, sm, maxdegree, small_k)
    return fval
    

def get_log_likelihood_window(frac1, frac2, window, k, pos, mu, pwm, rk):
    h = 1
    for j in rk:#range(k):
            h = h * pwm[window[j]][j]
    v = frac1 * np.exp(-math.pow((pos - mu), 2) / frac2) * h  
    
    return v


def compute_w(L, pwm, sigma, mu, k):
    frac1 = (1. / (math.sqrt(2 * math.pi) * sigma))
    frac2 = (2 * math.pow(sigma, 2))
    x = k - 1   
    po = int(math.pow(4, k))
    last_idx = k - 1
    w = np.zeros((po * L))
    window = np.zeros((k))
    rsig = int(round(sigma) + 0.5) #simga
    cival = 3 * rsig 
    ival = x + cival
    end = int(round(min(round(mu + ival), L - k + 1)))
    start = int(round(mu - ival))
    if mu < ival:
        start = 0
        
    if mu > L - ival:
        end = L - k + 1
    
    rk = range(k)
    r1 = range(int(start), int(end), 1)
    for i in range(po):
        for j in r1:#range(int(start), int(end), 1): 
            w[j * po + i] += get_log_likelihood_window(frac1, frac2, window, k, j, mu, pwm, rk)
              
        window[last_idx] += 1
        window_ptr = last_idx
        
        while window[window_ptr] == 4 and window_ptr > 0:
            window[window_ptr] = 0
            window_ptr -= 1
            window[window_ptr] += 1
            
    return w


def compute_scoring(mu, sigma, L, w, k, sm, D):  
    x = k - 1
    po = int(math.pow(4, k))
    div = 1. / po
    
    poim = np.zeros((po * L))
    
    v = []
    motifs = range(len(mu))
    for m in motifs:
        vhh = []
        for msub in range(D):
            vh = (np.zeros((po * (L + 2 * x))))
            vh[x * po:(L + x) * po] = w[m][msub]
            vhh.append(vh)
        v.append(vhh)
        
        
    zlcut = int(math.pow(2, k * 2) - 1)   
    
    M = len(mu)
    
    end = np.zeros((M, D))
    start = np.zeros((M, D))
    r2 = []
    for m in motifs:
        rsig = int(round(sigma[m]) + 0.5)
        cival = 3 * rsig
        ival = x + cival
        
        r2h = []
        for msub in range(D):
        
        
            end[m][msub] = int(round(mu[m] + msub + ival))
            start[m][msub] = int(round(mu[m] + msub - ival))
            if mu[m] + msub < ival:
                start[m][msub] = 0
            if mu[m] + msub > L - ival:
                end[m][msub] = L 
            r2h.append(range(int(start[m][msub]), int(end[m][msub]), 1))    
        r2.append(r2h)
            
              
    lenA = po * (2 * k - 1)
    AM = -np.ones(lenA)
    v1 = (k - 1) * po
    r1 = range(k - 1)
      
    
    for y in range(po):
        A = AM * div
        A[v1 + y] += 1       
        for i in r1:#range(k - 1):
            alloc = int(math.pow(4, i + 1))
            zlshift = int(math.pow(2, (k - i - 1) * 2))
            mov_bits = (i + 1) * 2
            zl = y << mov_bits & zlcut
            zr = y >> mov_bits
            
            comp_bits = math.pow(4, k - (i + 1))  
            '''
            h1: start with matrix allocation at position: [(k - 1) * po + y] which is the positional oligomer that
            we compare with. So the positional oligomers with (k-1) overlappings, which we have calculate with bit 
            wise shift, variable zl, to the right site lay in Matrix A at column k*po = (k+i)*po for the first i
            
            h2: is the columnvalue for the left side in A
            '''
            h1 = (k + i) * po + zl
            v2 = comp_bits * div
            A[h1:h1 + alloc] += v2
            h2 = (k - i - 2) * po + zr
            for pos in range(alloc): 
                A[ h2 + pos * zlshift] += v2
                
        for m in motifs:
            for msub in range(D):
                for j in r2[m][msub]:#range(int(start), int(end), 1): 
                    sp = j * po
                    poim[sp + y] += sm[m]*np.dot(v[m][msub][sp:(j + 2 * x + 1) * po], A)   
    return poim


def extend_func_POIM(Q, P, M, L, r, mu, sigma, sm, maxdegree, small_k):    
    
    R = []
    cR = []
    for k in range(small_k):
        cR.append(np.zeros((math.pow(4, k + 1) * L)))
        R.append(np.zeros((math.pow(4, k + 1), L)))
    
    
    for p in range(len(P)):
        v = []
        D = P[p] - small_k + 1
        for m in range(M[p]):
            vh = []
            for msub in range(D):
                vh.append(compute_w(L, r[p][m][:, msub:msub + small_k], sigma[p][m], mu[p][m] + msub, small_k))
            v.append(vh)
        cR[small_k - 1] += compute_scoring(mu[p], sigma[p], L, v, small_k, sm[p], D) #ein += statt = eingefuegt am 16.10
            
            
    for i in range(small_k):
        R[i] = utils.vec2matrix(cR[i], i + 1)  
    
    
    f = 0
    for p in range(len(P)):#small_k):
        for i in range(len(R[small_k - 1])): 
            for j in range(len(R[small_k - 1][0])):
                f += math.pow(R[small_k - 1][i][j] - Q[small_k - 1][i][j], 2) 
    fval = 1. / 2 * f     
    return fval
