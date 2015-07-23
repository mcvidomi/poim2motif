import copy
import math
import pdb
import random
import timeit
import cPickle as pickle
import numpy as np
from poim import *
import poim
from shogun.Features import *
from shogun.Kernel import *
from shogun.Classifier import *
from shogun.Evaluation import *
from shutil import *


dna = ['A', 'C', 'G', 'T']
def simulate_sequence(length):
    sequence = ''
    for i in range(length):
        sequence += random.choice(dna)
    return sequence


def mutate_motif(motiv,probmut):
    mutmot = ""
    for i in range(len(motiv)):
        rnd = random.random()
        if (rnd <= probmut):
            dnashort =['A', 'C', 'G', 'T']
            dnashort.pop(dnashort.index(motiv[i]))
            mutmot += random.choice(dnashort)
        else:
            mutmot +=motiv[i]
    return mutmot


def gensequences(tally,positives,sequenceno,prob,motif,mu):
    sequences = []
    labels = np.concatenate((ones(positives),-ones(sequenceno-positives)))
    #np.ones(sequenceno)*(-1)
    ml = len(motif)
    for i in range(sequenceno):
        aa = simulate_sequence(tally)
        if i < positives:
            #       labels[i]=1
            mut=mutate_motif(motif,prob)
            aa = aa.replace(aa[mu:mu + ml], mut)
        sequences.append(aa)


    return [sequences,labels]


def svmTraining(x,y,C,kernel_degree):
    feats_train  = StringCharFeatures(x,DNA)
    labels = BinaryLabels(y);
    start = timeit.default_timer()
    kernel = WeightedDegreePositionStringKernel(feats_train, feats_train, kernel_degree) 
    stop = timeit.default_timer()
    time_kernel = stop-start
    start = timeit.default_timer()
    svm = LibSVM(C, kernel, labels)    
    svm.train()
    stop=timeit.default_timer()
    time_svm = stop-start
    return [svm,time_kernel,time_svm]


#fm_train_dna = gensequences(tally,positives,sequenceno,0,motiv,mu)
def svmApply(svm,x):
    featstest = StringCharFeatures(x,DNA)
    outputs=svm.apply(featstest)
    #pm = PerformanceMeasures(labels_test, output);
    #acc = pm.get_accuracy();
    #roc = pm.get_auROC();
    #fms = pm.get_fmeasure();
    outputlabels = outputs.get_labels();
    return outputs.get_values(),outputlabels
 
def roc(y,scores,outputlabels,num_pos):
    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=num_pos)
    auc=metrics.auc(fpr, tpr)
    acc=metrics.accuracy_score(y,outputlabels)
    return fpr,tpr,auc,acc

    

def computePOIM(x,y,poim_degree,kernel_degree,savepath):
    feats_train  = StringCharFeatures(x,DNA)
    labels = BinaryLabels(y);
    kernel = WeightedDegreePositionStringKernel(feats_train, feats_train, kernel_degree)
    C=1
    svm = LibSVM(C, kernel, labels)    
    svm.train()
    tally=len(x[0])
    ma = poim.compute_poims(svm,kernel,poim_degree,tally)
    fobj = open(savepath,'wb')
    pickle.dump(ma,fobj)
    fobj.close()
    return ma  
