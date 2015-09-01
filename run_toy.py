import makeData as mkd
import makePOIM as mkp
import experiment
import pdb

# Task
TASK_1 = "compute POIM"
TASK_2 = "compute Motif"
TASK_3 = "plot POIM"

CURRENT_TASK = TASK_2

#read data
experiment_name = "ExpTOY1"
path = "/home/mvidovic/POIMall/experiments/"+ experiment_name + "/"
savepath = "/home/mvidovic/POIMall/data/toy/toy1.pkl"
poimpath = "/home/mvidovic/POIMall/data/toy/poim1.pkl"

if CURRENT_TASK == TASK_1:
    tally=30
    positives=250
    sequenceno=1000
    prob=0.0
    motif="GAGATT"
    mu=10
    x,y=mkd.gensequences(tally,positives,sequenceno,prob,motif,mu)
    #compute POIM
    poim_degree = 6
    kernel_degree = 10
    poims = mkp.computePOIM(x,y,poim_degree,kernel_degree,poimpath)


# Compute Motif
if CURRENT_TASK == TASK_2:
    
    if __name__ == '__main__':
        read_data = 1
        
        num_motifs=1
        maxdegree=6
        small_k=2
        
        optvar      = [1,1,1,0]
        ini_pwm     = "greedy"
        ini_mu      = [[10]]
        ini_sig     = [[1]]
        ini_lam     = [[1]]
        P           = [6]
        M           = [1]
        solver      = 'ralg'
        iprint      = 2
        maxIter     = 100
        ftol        = 1e-01
        gradtol     = 1e-01
        diffInt     = 1e-01
        contol      = 1e-01
        maxFunEvals = 1e03
        maxTime     = 1
        replacenucleotids = ""
        replaceposition = 0
        savepoims=0
        type = "real"
        exp = experiment.Experiment(type)
        
        exp.setup(path, poimpath, savepoims, optvar, ini_pwm, ini_mu, ini_sig, ini_lam, P, M,replacenucleotids,replaceposition, num_motifs, maxdegree, small_k, solver, iprint, maxIter, ftol, gradtol, diffInt, contol, maxFunEvals, maxTime)
        exp.parameter_prompt()
        exp.optimize()
