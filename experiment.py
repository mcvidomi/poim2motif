'''
Created on 08.06.2015

@author: marinavidovic
'''
import os
import utils
import motif
from datetime import datetime
import pdb
import numpy as np
from openopt import NLP
import func
import grad
import pickle
import view

class Experiment:
    '''
    classdocs
    '''
    num_motifs = 1
    poimdegree = 2
    lamall_k  = 2
    vecsize=np.zeros(4)
    optvar = np.ones(4)
    

    def __init__(self, type):
        '''
        type    ->     real or toy
        repetitions 
        '''
        self.type       = type 
        self.motif      = None
        self.poim       = None
        self.aQ         = None     
        self.Q          = None
        self.maxdegree  = None
        self.maxPOIM    = None
        self.diffPOIM   = None
        self.L          = None
        self.result     = None
        self.out_pwm    = None
        self.out_sig    = None
        self.out_lam    = None
        self.out_mu     = None
        self.out_fopt   = None
        self.out_time   = None
        self.out_iter   = None
        self.out_eavals = None
        self.out_istop  = None
        self.small_k   = None
        self.dictname   = None
        
        
    def setup(self,path,poimpath,savepoims,optvar,ini_pwm,ini_mu,ini_sig,ini_lam,P,M,replacenucleotids=[],replaceposition=[],num_motifs=1,maxdegree=6,small_k=2,solver='ralg',iprint=2, maxIter=10000, ftol=1e-03, gradtol=1e-03, diffInt=1e-03, contol=1e-03, maxFunEvals=1e04, maxTime=1000):
        self.num_motifs         = num_motifs
        self.small_k            = small_k
        self.maxdegree          = maxdegree
        self.path               = path
        self.poimpath           = poimpath
        self.motif              = motif.Motif()
        self.solver             = solver
        self.optvar             =optvar
        self.ini_mu             = ini_mu
        self.ini_sig            = ini_sig
        self.ini_lam             = ini_lam
        self.P                  = P 
        self.M                  = M 
        self.replacenucleotids  = replacenucleotids
        self.replaceposition    = replaceposition
        self.iprint             = iprint
        self.maxIter            = maxIter
        self.ftol               = ftol
        self.gradtol            = gradtol
        self.diffInt            = diffInt
        self.contol             = contol
        self.maxFunEvals        = maxFunEvals
        self.maxTime            = maxTime
        self.read(poimpath)
        self.iniPWM(ini_pwm)
        self.folderStruct(self.path)
        if savepoims:
                view.figurepoimsimple(self.Q[1], self.path + "/poim.eps", 0)
                view.figuremaxPOIM(self.diffPOIM, self.path + "/diffpoim.eps", 0)
    
    def iniPWM(self, ini_pwm):
        if ini_pwm == 'greedy':
            self.ini_pwm = utils.greedy_init(self.Q, self.P, self.M, self.ini_mu)
        elif ini_pwm == 'random':
            self.ini_pwm=utils.random_init(self.P, self.M)
        else:
            self.ini_pwm = ini_pwm 
            
    def read(self,poimpath):
            self.aQ, self.Q, self.maxdegree, self.maxPOIM, self.diffPOIM, self.L = utils.read_data_from_file(poimpath)
                  
    def folderStruct(self,path):
        dt = datetime.now()
        tt = dt.timetuple()
        self.dictname = "exp"
        for i in range(6):   
            self.dictname += "_"+str(tt[i])
        self.path = path+self.dictname
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        else:
            print "Experiment path already exists."    
          
          
    def input(self):
        '''
        organize optimization variable x regarding 
        optvar[0]: pwm
        optvar[1]: lam
        optvar[2]: sig
        optvar[3]: mu
        the values which are set to one in optvar are optimization variables
        '''
        x1, x2, x3, x4 = 0, 0, 0, 0
        x = []
        if (self.optvar[0] == 1):
            x1, self.vecsize[0] = utils.matr2vec(self.ini_pwm, self.P, self.M)
            x = np.append(x, x1)
        if (self.optvar[1] == 1):
            x2, self.vecsize[1] = utils.list2vec(self.ini_lam, self.P, self.M)
            self.vecsize[1] = self.vecsize[1] + self.vecsize[0]
        if (self.optvar[2] == 1):
            x3, self.vecsize[2] = utils.list2vec(self.ini_sig, self.P, self.M)
            self.vecsize[2] = self.vecsize[1] + self.vecsize[2]
        if (self.optvar[3] == 1):
            x4, self.vecsize[3] = utils.list2vec(self.ini_mu, self.P, self.M)
            self.vecsize[3] = self.vecsize[2] + self.vecsize[3]
        x = []
        if (self.vecsize[0] != 0):
            x = np.append(x, x1)
        if (self.vecsize[1] != 0):
            x = np.append(x, x2)
        if (self.vecsize[2] != 0):
            x = np.append(x, x3)
        if (self.vecsize[3] != 0):
            x = np.append(x, x4)
                
        lb_mu, lb_r_lam, lb_sig = [], [], []
        ub_mu, ub_r_lam, ub_sig = [], [], []
        lb_r_lam = np.ones(self.vecsize[1]) * 0.005
    
        if (self.vecsize[2] != 0):
            lb_sig = np.ones(self.vecsize[2] - self.vecsize[1]) * 0.0001
            ub_sig = np.ones(self.vecsize[2] - self.vecsize[1]) * (self.L - self.maxdegree + 1)
        if (self.vecsize[3] != 0):
            lb_mu = np.zeros(max(0, self.vecsize[3] - self.vecsize[2]))
            ub_mu = np.ones(max(0, self.vecsize[3] - self.vecsize[2])) * (self.L - self.maxdegree + 1)
        lb = np.append(np.append(lb_r_lam, lb_sig), lb_mu)
        ub_r_lam = np.ones(self.vecsize[1]) * 10     
        ub = np.append(np.append(ub_r_lam, ub_sig), ub_mu)
        
        
        lenA=int(np.max(self.vecsize))
        lenk=int(self.vecsize[0])/4
        Aeq=np.zeros((lenk,lenA))
        beq=np.ones(lenk)
        for i in range(lenk):
            for pk in range(i,lenA-2,lenk):
                Aeq[i,pk] = 1
        return x,Aeq,beq,lb,ub
          
    def optimize(self):
        x,Aeq,beq,lb,ub=self.input()
        #p = NLP(func.extend_callfunc, x, df=grad.extend_callgrad, Aeq=Aeq,beq=beq, lb=lb, ub=ub, args=(self.Q, self.P, self.M, self.L,self.ini_pwm, self.ini_mu, self.ini_sig, self.ini_lam, self.maxdegree, self.vecsize, "5", self.optvar, self.lamall_k),  diffInt=self.diffInt, ftol=self.ftol, plot=0, iprint=self.iprint,maxIter = self.maxIter, maxFunEvals = self.maxFunEvals, show=False, contol=self.contol)
        self.parameter_prompt()
        
        p = NLP(func.extend_callfunc, x, df=grad.extend_callgrad, lb=lb, ub=ub, args=(self.Q, self.P, self.M, self.L,self.ini_pwm, self.ini_mu, self.ini_sig, self.ini_lam, self.maxdegree, self.vecsize, "5", self.optvar, self.small_k),  diffInt=self.diffInt, ftol=self.ftol, plot=0, iprint=self.iprint,maxIter = self.maxIter, maxFunEvals = self.maxFunEvals, show=False, contol=self.contol)
        #p.checkdf()
        #p.checkdc()
        #p.checkdh()
        self.result = p._solve(self.solver)
        self.output()
        
   
    
    def output(self):
        
        #==========================================================================
        # 
        #p.iterTime = []
        # p.iterValues.x = [] # iter points
        # p.iterValues.f = [] # iter ObjFunc Values
        # p.iterValues.r = [] # iter MaxResidual
        # p.iterValues.rt = [] # iter MaxResidual Type: 'c', 'h', 'lb' etc
        # p.iterValues.ri = [] # iter MaxResidual Index
        # p.solutions = [] # list of solutions, may contain several elements for interalg and mb other solvers
        #==========================================================================
        
        
        '''
        save raw result vector
        '''
        fobj = open(self.path+ "/result_raw_" + self.dictname + ".pkl", 'wb')
        pickle.dump(self.result, fobj)
        fobj.close()
        
        
        if (self.optvar[0] == 1):
            self.out_pwm = utils.converter(self.result.xf[:self.vecsize[0]], self.P, self.M)
        else:
            self.out_pwm = self.ini_pwm
        if (self.optvar[1] == 1):
            self.out_lam = utils.listconverter(self.result.xf[self.vecsize[0] :self.vecsize[1]], self.P, self.M)
        else:
            self.out_lam = self.ini_lam
        if (self.optvar[2] == 1):
            self.out_sig = utils.listconverter(self.result.xf[self.vecsize[1]:self.vecsize[2]], self.P, self.M)
        else:
            self.out_sig = self.ini_sig
        if (self.optvar[3] == 1):
            self.out_mu = utils.listconverter(self.result.xf[self.vecsize[2]:], self.P, self.M)
        else:
            self.out_mu = self.ini_mu
      
        self.out_fopt   = self.result.ff
        self.out_time   = self.result.elapsed['solver_time']
        self.out_iter   = None
        self.out_eavals = self.result.evals['f']
        self.out_istop  = self.result.istop
        
        fobj = open(self.path+ "/output_"+self.dictname+".pkl", 'wb')
        pickle.dump([self.out_pwm, self.out_lam, self.out_sig, self.out_mu,self.out_fopt, self.out_time,self.out_eavals,self.out_iter,self.out_istop ], fobj)
        fobj.close()
        
        utils.parameter2file(self.path + "/result_prepared_" + self.dictname + ".txt",self.solver,self.small_k,self.ini_pwm,self.ini_lam,self.ini_sig,self.ini_mu,self.out_pwm, self.out_lam, self.out_sig, self.out_mu,self.out_fopt, self.out_time,self.out_eavals,self.out_iter,self.out_istop,self.optvar, self.M, self.P,self.maxIter,self.ftol,self.gradtol,self.diffInt,self.contol,self.maxFunEvals,self.maxTime)
        utils.makeJasparMotif(self.out_pwm, 1, 1, self.P, self.M, self.path+ "/jaspar_"+self.dictname+".txt",self.replacenucleotids,self.replaceposition)
        utils.seqLogo(self.path,"/jaspar_"+self.dictname)
        
        
    def parameter_prompt(self):
        utils.parameter_prompt(self.Q, self.ini_pwm, self.ini_mu, self.ini_sig, self.ini_lam, self.M, self.P, self.L, self.diffPOIM, "optimization parameters")
        
            
