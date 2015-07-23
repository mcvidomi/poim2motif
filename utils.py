'''
Created on 08.06.2015

@author: marinavidovic
'''

import pickle 
import math
import numpy as np
import random
import pdb
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

def seqLogo(path2pwm,filename):
    path=path2pwm + filename +".txt"
    importr("grid")
    importr("seqLogo")
    r = robjects.r
    m = r('''m = read.table( ( "%(path)s"))''' % locals() )
    pwm=r.makePWM(m)
    s = r.paste(path2pwm,filename,"_seqLogo.pdf",sep="")
    r.pdf(file=s)
    r.seqLogo(pwm,r("ic=FALSE"))
    r("dev.off()")

def dna2vec(dna):
    """
    converts the digits to the dna alphabet and returns dna string
    """
    str = []
    for i in range(len(dna)):
        if dna[i] == 'A':
            str.append(0)
        elif dna[i] == 'C':
            str.append(1)
        elif dna[i] == 'G':
            str.append(2)
        else:
            str.append(3)
    return str

def loopidx2dna(index, motivelen):
    dna = []
    dnah = []
    for j in range(motivelen):
        dnah.append(index % 4)
        index = int(index / 4)
    dnah.reverse()
    dna.append(vectodna(dnah))  
    
    return dna, dnah

def vectodna(vec):
    """
    converts the digits to the dna alphabet and returns dna string
    """
    str = []
    for i in range(len(vec)):
        if vec[i] == 0:
            str.append('A')
        elif vec[i] == 1:
            str.append('C')
        elif vec[i] == 2:
            str.append('G')
        else:
            str.append('T')
    return str

def veclisttodna(vec):
    """
    converts the digits to the dna alphabet and returns dna string
    """
    str = ""
    for i in range(len(vec)):
        if vec[i] == 0:
            str += 'A'
        elif vec[i] == 1:
            str += ('C')
        elif vec[i] == 2:
            str += ('G')
        else:
            str += ('T')
    return str

def vec2matrix(vec, k):
    po = math.pow(4, k)
    tally = int(len(vec) / po)
    matr = np.zeros((po, tally))
    for l in range(tally):
            matr[:, l] = vec[l * po:l * po + po]
    return matr

def motif_extraction(pwm, mu, P, M):
    motif = []
    dna = ['A', 'C', 'G', 'T']
    
    for p in range(len(P)):
        motifm = []
        for m in range(M[p]):
            idx = []
            for k in range(len(pwm[p][m][0])):
                maxval = max(pwm[p][m, :, k])
                idx.append(dna[(pwm[p][m, :, k]).tolist().index(maxval)])
            
            motifm.append(idx)  
        motif.append(motifm)
    return motif

def converter(par, P, M):
    if type(par[0]) == np.float64 :
        r = []
        count = 0
        for p in range(len(P)):
            pwm = np.zeros((M[p], 4, P[p]))
            for m in range(M[p]):
                for i in range(4):
                    for j in range(P[p]):
                        pwm[m][i][j] = par[count]
                        count += 1
                        
            r.append(pwm)
        return r
    return par

def listconverter(par, P, M):
    if type(par[0]) == np.float64 :
        vec = []
        count = 0
        for p in range(len(P)):
            vecm = []
            for m in range(M[p]):
                vecm.append(par[count])
                count += 1
            vec.append(vecm)
        return vec
        
    return par


def matr2vec(matr, P, M):
    add = 0
    for p in range(len(P)):
        add += P[p]*M[p]*4
    vec = np.zeros(add)
    count = 0
    for p in range(len(P)):
        for m in range(M[p]):
            for i in range(4):
                vec[count:count + P[p]] = matr[p][m][i]
                count += P[p]
    return vec, add


def list2vec(list, P, M):
    add = 0
    for p in range(len(P)):
        add += M[p]
    vec = np.zeros(add)
    count = 0
    for p in range(len(P)):
        vec[count:count + M[p]] = list[p]
        count += M[p]
    return vec, add

def round_pwm(pwm, P, M):
    r = []
    
    for p in range(len(P)):
        pwmh = np.zeros((M[p], 4, P[p]))
        for m in range(M[p]):
            for i in range(4):
                for j in range(len(pwm[p][m][0])):
                    pwmh[m][i][j] = round(pwm[p][m][i][j], 2)
                                        
        r.append(pwmh)
    return r

def read_data_from_file(POIMfile):

    fobj = open(POIMfile, 'rb')
    
    aQ = pickle.load(fobj)
    Q = aQ[0]
    maxdegree = len(Q)
    maxPOIM = aQ[1]
    diffPOIM = aQ[2]
    L = len(Q[1][0])
    return aQ, Q, maxdegree, maxPOIM, diffPOIM, L

  
def startpwmrnd(motivelen, motivno):
    '''
    generates a random start pwm for a motive with length motivelength
    '''
    pwm = np.zeros((motivno, 4, motivelen))
    for m in range(motivno):
        for j in range(int(motivelen)):  
            sumcol = 0
            for i in range(4):
                pwm[m][i][j] = random.uniform(1, 100)
                sumcol += pwm[m][i][j]
            for i in range(4):
                pwm[m][i][j] = round(pwm[m][i][j] / sumcol, 2)
    return pwm


def best_pwm_ini(motivno, motivelen, exp_motif):   
    pwm = np.ones((motivno, 4, motivelen)) * 0.1
    for m in range(motivno):
        for j in range(int(motivelen)):  
            pwm[m][exp_motif[m][j]][j] = 0.7
    return pwm

def greedy_init(Q, P, M, mu):
    r = [] 
    mo, mop = POIM22motif(Q, P, M, mu)
    for d in range(len(P)):
        r.append(best_pwm_ini(M[d], P[d], mop[d]))
    return r

def random_init(P, M):
    
    r = []
    for d in range(len(P)):
        r.append(startpwmrnd(P[d], M[d]))
    return r

def POIM22motif(Q, P, M, givenmu):
    mo, mop = [],[]

    for p in range(len(P)):
        mohh, moph = [],[]
        for m in range(M[p]):
            moh = ""
            for e in range(P[p]):
                maxval = max(Q[1][:, givenmu[p][m] + e])
                idx = (Q[1][:, givenmu[p][m] + e]).tolist().index(maxval)
                rm, rr = loopidx2dna(idx, 2)
                moh = moh + rm[0][0]       
            mohh.append(moh)
            moph.append(dna2vec(moh))
        mop.append(moph)
        mo.append(mohh)
    return mo, mop

def parameter_prompt(Q, r, mu, sigma, lam, M, P, L, maxPOIM, task_prompt):
    print "\n\n--------------------------------------------------------------------------------------"
    print  task_prompt , ":\n"
    print '    {0:10} ==> {1:10}\n'.format('motif length ', P)
    print '    {0:10} ==> {1:10}\n'.format('number of motifs ', M)
    print '    {0:10} ==> {1:10}\n'.format('sequence length:   ', L)

    motifs = motif_extraction(r, mu, P, M)
    
    r = round_pwm(r, P, M)    
    
    for p in range(len(P)):
        for m in range(M[p]):
            print "\n\n     ", m + 1, ". motif of length :" , P[p], "\n"
            print '         lam =     ', lam[p][m]
            print '         sigma =  ', sigma[p][m]
            print '         mu =     ', mu[p][m]
            print '         pwm =    '
            
            for i in range(4):
                print "                  ", r[p][m][i]    
            print '         motif =  ', motifs[p][m]

    print "--------------------------------------------------------------------------------------\n\n"
    
def parameter2file(path,solver,lamall_k,ini_pwm, ini_lam, ini_sig, ini_mu,out_pwm, out_lam, out_sig, out_mu,out_fopt, out_time,out_eavals,out_iter,out_istop,optvar, M, P,maxIter,ftol,gradtol,diffInt,contol,maxFunEvals,maxTime):
    fobj = open(path, 'w')
    fobj.write("optimization parameter [pwm,lambda,sigma,mu]: \n" + str(optvar))
    fobj.write("solver = " + str(solver) + "\n") 
    fobj.write("maxIter = " + str(maxIter) + "\n")
    fobj.write("ftol = " + str(ftol) + "\n")
    fobj.write("gradtol = " + str(gradtol) + "\n")
    fobj.write("diffInt = " + str(diffInt) + "\n")
    fobj.write("contol = " + str(contol) + "\n")
    fobj.write("maxFunEvals = " + str(maxFunEvals) + "\n")
    fobj.write("maxTime = " + str(maxTime))
    fobj.write("\n")
    fobj.write("lamall k = " + str(lamall_k))
    fobj.write("\n")    
      
    dna = ['A', 'C', 'G', 'T']
    for p in range(len(P)):
        fobj.write("\n\n")
        for m in range(M[p]):
            fobj.write("Motif " + str(m + 1) + "    length = " + str(P[p]) + "\n")
            
        
            fobj.write('input parameters:  \n\n')
            fobj.write('pwm: ' + str(ini_pwm[p][m]) + '\n')
            fobj.write('lambda: ' + str(ini_lam[p][m]) + '\n')
            fobj.write('sigma: ' + str(ini_sig[p][m]) + '\n')
            fobj.write('mu: ' + str(ini_mu[p][m]) + '\n')
            
            fobj.write('\n\noutput parameters:  \n\n')
            fobj.write('fopt: '   + str(out_fopt) + '\n')
            fobj.write('time: '    + str(out_time) + '\n')
            fobj.write('fevals: '       + str(out_eavals) + '\n')
            fobj.write('fiter: '   + str(out_iter) + '\n')
            fobj.write('istop: '   + str(out_istop) + '\n')
            fobj.write('lambda: '    + str(out_lam[p][m]) + '\n')
            fobj.write('sigma: '    + str(out_sig[p][m]) + '\n')
            fobj.write('mu: '       + str(out_mu[p][m]) + '\n')
            fobj.write("\n")
             
            #pwm output
            for row in range(4):
                fobj.write(dna[row] + ' ')
                
                fobj.write(str(out_pwm[p][m][row])+'\n')
              
            fobj.write('\n\n')

def makeJasparMotif(out_pwm, pot, normon, P, M, outputJaspar,replacenucleotids,replaceposition):
    pwm=pwmnorm(out_pwm[0][0], normon)
    #insert non-polymorphic loci 
    for i in range(len(replacenucleotids)):
        pwm[:,replaceposition+i]=0
        if replacenucleotids[i]=='A':
            pwm[0,replaceposition+i]=1
        elif replacenucleotids[i]=='C':
            pwm[1,replaceposition+i]=1
        elif replacenucleotids[i]=='G':
            pwm[2,replaceposition+i]=1        
        else:
            pwm[3,replaceposition+i]=1            
    np.savetxt(outputJaspar, pwm)
    
def pwmnorm(r,norm_on):
    """
    with the little hack:
    after normalization there are too many comma values - but round them of about 2 comma values, 
    the column total is not always 1, so the function adds the difference from 1 and the column total to 
    A at the column.
    originally it was: r_normh[m][j][i] = rh[j][i]*div
    """ 
    l =len(r[0])
    r_norm = np.zeros((4,l))
    r_sum = sum(r)
    for i in range(len(r[0])):
        div = (norm_on +0.0) / r_sum[i]
        for j in range(4):
            r_norm[j][i] = r[j][i]*div        
    return r_norm
