'''
Created on 02.07.2012

@author: Marina
'''



import math
import utils
import numpy as np
import copy
import pickle
import pdb


def extend_callgrad(x, Q, P, M, L, r, mu, sigma, sm, maxdegree, vecsize, experiment, opt_var, small_k):
    #opt_var elemnts
    #         1. pwm
    #         2. sm
    #         3. sigma
    # 0: no optimization variable
    # 1: optimization variable
    #example: opt_var =[1,0,1] === pwm,sigma
    
    if (opt_var[0] == 1):
        r = utils.converter(x[:vecsize[0]], P, M)
    if (opt_var[1] == 1):
        sm = utils.listconverter(x[vecsize[0] :vecsize[1]], P, M)
    if (opt_var[2] == 1):
        sigma = utils.listconverter(x[vecsize[1]:vecsize[2]], P, M)
    if (opt_var[3] == 1):
        mu = utils.listconverter(x[vecsize[2]:], P, M)
    
    gradmu, gradsig, gradpwm, gradsm = extend_gradient(Q, P, M, L, r, mu, sigma, sm, maxdegree, opt_var, small_k)
    
    gradient = []
    if (opt_var[0] == 1):
        gradpwm, vecsize1 = utils.matr2vec(gradpwm, P, M)
        gradient = np.concatenate((gradient, gradpwm))
    if (opt_var[1] == 1):
        gradsm, vecsize2 = utils.list2vec(gradsm, P, M)
        gradient = np.concatenate((gradient, gradsm))
    if (opt_var[2] == 1):
        gradsig, vecsize3 = utils.list2vec(gradsig, P, M)
        gradient = np.concatenate((gradient, gradsig))
    if (opt_var[3] == 1):
        gradmu, vecsize4 = utils.list2vec(gradmu, P, M)
        gradient = np.concatenate((gradient, gradmu))
    return gradient


def dh(x, args):
    return 2 * (x - 5)

def get_log_likelihood_window(frac1, frac2, window, k, pos, mu, pwm, wgradpwm, t, rk):
    h = 1
    for j in rk:
            h = h * pwm[window[j]][j]
            wgradpwm[window[j]][j][t] = 1. / pwm[window[j]][j]
    v = frac1 * np.exp(-math.pow((pos - mu), 2) / frac2) * h  
    
    return v


def compute_w(L, pwm, sigma, mu, k):
    frac1 = (1. / (math.sqrt(2 * math.pi) * sigma))
    frac2 = (2 * math.pow(sigma, 2))
    frac3 = math.pow(sigma, 2)
    frac4 = math.pow(sigma, 3)
    x = k - 1   
    po = int(math.pow(4, k))
    last_idx = k - 1
    w = np.zeros((po * L))
    wgradmu = np.zeros((po * L))
    wgradsig = np.zeros((po * L))
    wgradpwm = np.zeros((4, k, po * L))
    window = np.zeros((k))
    rsig = int(round(sigma) + 0.5) #sigma#
    cival = 3 * rsig 
    ival = x + cival
    end = int(round(min(mu + ival, L - k + 1)))
    start = int(round(mu - ival))
    
    if mu < ival:
        start = 0
        
    if mu > L - ival:
        end = L - k + 1
        
    rk = range(k)
    r1 = range(int(start), int(end), 1)
    

    for i in range(po):            
        for j in r1:
            t = j * po + i            
            w[t] += get_log_likelihood_window(frac1, frac2, window, k, j, mu, pwm, wgradpwm, t, rk)
            wgradmu[t] += copy.deepcopy(w[t]) * (j - mu) / frac3
            wgradsig[t] += copy.deepcopy(w[t]) * (math.pow(j - mu, 2) / frac4 - 1. / sigma)
            
            
            
        window[last_idx] += 1
        window_ptr = last_idx
        
        while window[window_ptr] == 4 and window_ptr > 0:
            window[window_ptr] = 0
            window_ptr -= 1
            window[window_ptr] += 1        
    return w, wgradmu, wgradsig, wgradpwm



def compute_scoring(P, M, L, r, mu, sigma, maxdegree, w, wgradmu, wgradsig, wgradpwm, k, sm, opt_var, D):
    x = k - 1
    po = int(math.pow(4, k))
    div = 1. / po
    

    K = D + k - 1

    R = []
    Rgradmu = []
    Rgradsig = []
    Rgradpwm = []
    v = []
    vgradmu = []
    vgradsig = []
    vgradpwm = []
    
    motifs = range(len(mu))
    for m in motifs:
        vhh = []
        vgradmuhh = []
        vgradsighh = []
        vgradpwmhh = []
        
        for msub in range(D):
            
            
            R.append(np.zeros((po * L)))
            Rgradmu.append(np.zeros((po * L)))
            if opt_var[2] != 0:
                Rgradsig.append(np.zeros((po * L)))
            Rgradpwm.append(np.zeros((4, K, po * L)))
            vh = np.zeros((po * (L + 2 * x)))
            vgradmuh = np.zeros((po * (L + 2 * x)))
            vgradsigh = np.zeros((po * (L + 2 * x)))
            vgradpwmh = np.zeros((4, K, po * (L + 2 * x)))
            vh[x * po:(L + x) * po] = w[m][msub]
            
            vgradmuh[x * po:(L + x) * po] = wgradmu[m][msub]
            vgradsigh[x * po:(L + x) * po] = wgradsig[m][msub]
            
            vhh.append(vh)
            vgradmuhh.append(vgradmuh)
            vgradsighh.append(vgradsigh)
            
            for n in range(4):
                for u in range(k):
                    vgradpwmh[n][u + msub][x * po:(L + x) * po] = wgradpwm[m][msub][n][u]
            vgradpwmhh.append(vgradpwmh)
            
            
        v.append(vhh)

        vgradmu.append(vgradmuhh)
        vgradsig.append(vgradsighh)
        vgradpwm.append(vgradpwmhh)
        
        
 
        
        
    zlcut = int(math.pow(2, k * 2) - 1)   
 
    end = np.zeros((len(mu), D))
    start = np.zeros((len(mu), D))
   
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
        
        
    AM = -np.ones(po * (2 * k - 1))
    v1 = (k - 1) * po
   
        

    for y in range(po):
        A = AM * div
        A[v1 + y] += 1       
        for i in range(k - 1):
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
            r5 = range(alloc)
            for pos in r5:
                A[ h2 + pos * zlshift] += v2
        
        for m in motifs:
            
            for msub in range(D):
            
                for j in r2[m][msub]:
                      
                    sp = j * po
                    R[m][sp + y] += np.dot(v[m][msub][sp:(j + 2 * x + 1) * po], A)      #sm[m]*
                    #Rgradmu[m][sp + y] += np.dot(vgradmu[m][sp:(j + 2 * x + 1) * po], A)  #
                    
                    Rgradmu[m][sp + y] += np.dot(vgradmu[m][msub][sp:(j + 2 * x + 1) * po], A)
                    
                    if opt_var[2] != 0:
                        Rgradsig[m][sp + y] += np.dot(vgradsig[m][msub][sp:(j + 2 * x + 1) * po], A)  
                    for n in range(4):
                        for u in range(k):
                            Rgradpwm[m][n][u + msub][sp + y] += np.dot(v[m][msub][sp:(j + 2 * x + 1) * po]*vgradpwm[m][msub][n][u + msub][sp:(j + 2 * x + 1) * po], A) 

    return R, Rgradmu, Rgradsig, Rgradpwm



def extend_gradient(Q, P, M, L, r, mu, sigma, sm, maxdegree, opt_var, small_k):
    R = []
    cR = []
    for k in range(small_k):
        cR.append(np.zeros((math.pow(4, k + 1) * L)))
        R.append(np.zeros((math.pow(4, k + 1), L)))
    
    
    
    
    
    Rgradmu = []  
    Rgradsig = []  
    Rgradpwm = []
    Rgradsm = []
    for p in range(len(P)): 
        D = P[p] - small_k + 1
        K = P[p]
        Rgradm = []
        if opt_var[2] != 0:
            Rgradsigm = []
        Rgradpwmm = []
        Rgradsmm = []
        
        v = []
        vgrad = []
        vgradsig = []
        vgradpwm = []
        for m in range(M[p]):
            vhh = []
            vgradhh = []
            vgradsighh = []
            vgradpwmhh = []
            
            for msub in range(D):
            
                vh, vgradh, vgradsigh, vgradpwmh = compute_w(L, r[p][m][:, msub:msub + small_k], sigma[p][m], mu[p][m] + msub, small_k)
                vhh.append(vh)
                vgradhh.append(vgradh)
                vgradsighh.append(vgradsigh)
                vgradpwmhh.append(vgradpwmh)
                
            v.append(vhh)
            vgrad.append(vgradhh)
            vgradsig.append(vgradsighh)
            vgradpwm.append(vgradpwmhh)
        
            
        R1, R2, R3, R4 = compute_scoring(P, M, L, r[p], mu[p], sigma[p], maxdegree, v, vgrad, vgradsig, vgradpwm, small_k, sm[p], opt_var, D)

        for m in range(M[p]):
            Rgradm.append(utils.vec2matrix(sm[p][m]*R2[m], small_k))
            if opt_var[2] != 0:
                Rgradsigm.append(utils.vec2matrix(sm[p][m]*R3[m], small_k))
            if opt_var[1] != 0:
                Rgradsmm.append(utils.vec2matrix(R1[m], small_k))
            for n in range(4):
                for u in range(P[p]):  
                    Rgradpwmm.append(utils.vec2matrix(sm[p][m]*R4[m][n][u], small_k))
            cR[small_k - 1] += sm[p][m]*R1[m]
            

        Rgradmu.append(Rgradm)
        if opt_var[2] != 0:
            Rgradsig.append(Rgradsigm)
        Rgradpwm.append(Rgradpwmm)
        if opt_var[1] != 0:
            Rgradsm.append(Rgradsmm)
    
    for i in range(small_k):
        R[i] = utils.vec2matrix(cR[i], i + 1)  

    
    
    
    gradmu = []
    gradsig = []
    gradpwm = []
    gradsm = []
    
    for d in range(len(P)):
        gradmu.append([0.]* M[d])
        gradsig.append([0.]* M[d])
        gradpwm.append(np.zeros((M[d], 4, P[d])))
        gradsm.append([0.]* M[d])
       
    ### here is the inner derivative multiplied with the outer derivative    
    for p in range(len(P)):
        for i in range(len(R[small_k - 1])): 
            for j in range(len(R[small_k - 1][0])):
                diff = R[small_k - 1][i][j] - Q[small_k - 1][i][j]
                count = 0
                for m in range(M[p]):
                    gradmu[p][m] += (diff * Rgradmu[p][m][i][j])
                    if opt_var[2] != 0:
                        gradsig[p][m] += (diff * Rgradsig[p][m][i][j])
                    if opt_var[1] != 0:
                        gradsm[p][m] += (diff * Rgradsm[p][m][i][j])
                
                    for n in range(4):
                        for u in range(P[p]): 
                            gradpwm[p][m][n][u] += (diff * Rgradpwm[p][count][i][j])
                            count += 1
    fobj2 = open("pruef_grad1", "wb")
    pickle.dump([gradmu, gradsig, gradpwm, gradsm], fobj2)
    fobj2.close()
    return gradmu, gradsig, gradpwm, gradsm
