'''
    Created on 08.06.2015
    
    @author: marinavidovic
    '''
import os
import pdb
import utils_svm
import pickle
import numpy as np
import copy
import genQ
import makePOIM
import view
import matplotlib
matplotlib.use('Agg')


if __name__ == '__main__':
    
    read_data = 1
    datapath = "/home/mvidovic/POIMall/data/real/human_acceptor_splice_data.txt"
    savepath = "/home/mvidovic/POIMall/data/real/human_acceptor_splice_data0.pkl"

    lines=1000
    if read_data:
        x,y=utils_svm.extractRealData(datapath,savepath,lines)
    else:
        fobj=open(savepath,'rb')
        x,y=pickle.load(fobj)
        fobj.close()
    num_pos = 100
    num_neg = 4*num_pos
    print "reduce samples"
    x_red,y_red = utils_svm.reduce_samples(copy.deepcopy(x),copy.deepcopy(y),num_pos,num_neg)
    nploci_letters,nploci_positions = utils_svm.non_polymorphic_loci(x_red)
    #read data
    experiment_name = "real1"
    if not os.path.exists(experiment_name):
        os.makedirs(experiment_name)
    poimpath=experiment_name+"/poim.pkl"
    tally=30
    positives=25
    sequenceno=100
    mutation_prob=0.0
    motif="ATTTT"
    mu=13
    x,y = makePOIM.gensequences(tally,positives,sequenceno,mutation_prob,motif,mu)
    
    #compute POIM
    poim_degree = 6
    kernel_degree = 8
    print "start poim computation"
    poims = makePOIM.computePOIM(x,y,poim_degree,kernel_degree,poimpath)
    Q2 = poims[0][1]
    #view.test()
    view.figurepoimsimple(Q2, "poim_pic", 0)
