import makeData as mkd
import makePOIM as mkp
import experiment
import pdb


# Task
TASK_1 = "compute POIM"
TASK_2 = "compute Motif"
TASK_3 = "plot POIM"

CURRENT_TASK = TASK_2


# Read Data
experiment_name = "ExpSigma"
path = "/home/mvidovic/POIMall/experiments/"+ experiment_name + "/"
datapath = "/home/mvidovic/POIMall/data/real/human_acceptor_splice_data0.txt"
savepath = "/home/mvidovic/POIMall/data/real/human_acceptor_splice_data0.pkl"
lines=1000
poimpath = "/home/mvidovic/POIMall/data/real/human_acceptor_splice_data0_poim_lines" + str(lines) +  ".pkl"

# Compute POIM
if CURRENT_TASK == TASK_1:
    x,y=mkd.extractRealData(datapath,savepath,lines)
    nploci = mkd.non_polymorphic_loci(x)
    
    poim_degree = 6
    kernel_degree = 10
    poims = mkp.computePOIM(x,y,poim_degree,kernel_degree,poimpath)

# Compute motif
if CURRENT_TASK == TASK_2:
    if __name__ == '__main__':

        read_data = 1
        num_motifs=1
        maxdegree=6
        small_k=2
        
        optvar      = [1,1,1,0]
        ini_pwm     = "greedy"
        ini_mu      = [[45]]
        ini_sig     = [[1]]
        ini_lam     = [[1]]
        P           = [20]
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
        replacenucleotids = "AG"
        replaceposition = 15
        savepoims=0
        type = "real"
        exp = experiment.Experiment(type)
        
        exp.setup(path, poimpath, savepoims, optvar, ini_pwm, ini_mu, ini_sig, ini_lam, P, M,replacenucleotids,replaceposition, num_motifs, maxdegree, small_k, solver, iprint, maxIter, ftol, gradtol, diffInt, contol, maxFunEvals, maxTime)
        exp.parameter_prompt()
        exp.optimize()




