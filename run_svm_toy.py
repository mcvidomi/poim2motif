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
    
    
    #datapath = "/Users/marinavidovic/Documents/work/POIM/data/human_acceptor_splice_data.txt"
    #savepath = "/Users/marinavidovic/Documents/work/POIM/data/human_acceptor_splice_data_xy.pkl"
    #poimfolder = "/Users/marinavidovic/Documents/work/POIM/data/human_acceptor_splice_data"
    #poimpath = "/Users/marinavidovic/Documents/work/POIM/data/human_acceptor_splice_data"+str(lines)+".pkl"
    #path = "/Users/marinavidovic/Documents/work/POIM/experiments/" + experiment_name + "/"
    
    #read data
    experiment_name = "toy1"
    if not os.path.exists(experiment_name):
        os.makedirs(experiment_name)
    poimpath=experiment_name+"poim.pkl"
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
