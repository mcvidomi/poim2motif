'''
    Created on 08.06.2015
    
    @author: marinavidovic
    '''
import os
import experiment
import pdb
import view
import utils

if __name__ == '__main__':
    
    
    experiment_name="real1"
    
    #read data
    #datapath = "/home/mvidovic/POIMall/data/real/human_acceptor_splice_data0.txt"
    
    
    poimpath = experiment_name + "/poim.pkl"
    path = experiment_name + "/"
    num_motifs=1
    maxdegree=6
    lamall_k=2
    
    optvar      = [1,1,1,0]
    ini_pwm     = "greedy"
    ini_mu      = [[10]]
    ini_sig     = [[1]]
    ini_lam      = [[1]]
    P           = [6]
    M           = [1]
    solver      = 'ralg'
    iprint      = 2
    maxIter     = 100
    ftol        = 1e-03
    gradtol     = 1e-01
    diffInt     = 1e-01
    contol      = 1e-01
    maxFunEvals = 1e03
    maxTime     = 1
    savepoims=1
    replacenucleotids=""
    replaceposition=0
    type = "real"
    exp = experiment.Experiment(type)
    exp.setup(path, poimpath, savepoims, optvar, ini_pwm, ini_mu, ini_sig, ini_lam, P, M,replacenucleotids,replaceposition, num_motifs, maxdegree, lamall_k, solver, iprint, maxIter, ftol, gradtol, diffInt, contol, maxFunEvals, maxTime)
    exp.parameter_prompt()
    
    pdb.set_trace()
    
    exp.optimize()
    
    pdb.set_trace()