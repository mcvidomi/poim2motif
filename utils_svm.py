import random
import pickle
import numpy as np
import copy
import math
import pdb
from shutil import *
import timeit
#from sklearn import metrics
#import sklearn
import poim
from shogun.Classifier import *
from shogun.Evaluation import *
from shogun.Features import *
from shogun.Kernel import *


def simulate_sequence(length):
    sequence = ''
    dna = ['A', 'C', 'G', 'T']
    for i in range(length):
        sequence += random.choice(dna)  #zufaellig Element aus dna anhaengen
    return sequence


def mutate_motif(motiv,probmut):
    dna = ['A', 'C', 'G', 'T']
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
    ml = len(motif)
    for i in range(sequenceno):
        aa = simulate_sequence(tally)
        if i < positives:
            
            mut=mutate_motif(motif,prob)
            aa = aa.replace(aa[mu:mu + ml], mut)
        sequences.append(aa)
    return sequences


def non_polymorphic_loci(x):
    counter = np.zeros((4,len(x[0])))
    for i in range(len(x)):
        for j in range(len(x[0])):
            if x[i][j] ==  'A':
                    counter[0,j]=counter[0,j]+1
            elif x[i][j] ==  'C':  
                    counter[1,j]=counter[1,j]+1
            elif x[i][j] ==  'G':  
                    counter[2,j]=counter[2,j]+1
            else:  
                counter[3,j]=counter[3,j]+1           
        
    counter=counter/len(x)
    
    dna = ['A', 'C', 'G', 'T']
    letter = []
    position = []
    for i in range(len(counter[0])):
        for j in range(4):
            if counter[j,i] == 1.0:
                print "nucleotid " , dna[j]," position", str(i)  
                letter.append(dna[j])
                position.append(i)
    

    return letter,position


def extractRealData(datapath,savepath,lines):
    data =file(datapath).readlines()[:lines]
    labels = []
    x = []
    cn=0
    for i in range(len(data)):
        labels.append(int(data[i][0:2]))
        x.append(data[i][3:-1])
        if int(data[i][0:2])==1:
            cn=cn+1
    print "numper of positive labels: " , cn
    if savepath !="":
        x=np.array(x)
        labels=np.array(labels)
        fobj = open(savepath,'wb')
        pickle.dump([x,labels],fobj) 
        fobj.close()
    return x,labels

def reduce_samples(x,y,num_pos,num_neg):
    xpos = x[y== 1]
    xneg = x[y==-1]
    if num_pos>len(xpos):
        print "Number of positive samples " + str(len(xpos)) + " is smaller than num_pos " +str(num_pos) +". Set num_pos to maximal available positive samples"
        num_pos = len(xpos)
    if num_neg>len(xneg):
        print "Number of negative samples " + str(len(xpos)) + " is smaller than num_neg " +str(num_pos) +". Set num_neg to maximal available negative samples"
        num_neg = len(xneg)
    x_red = xpos[0:num_pos].tolist() + xneg[0:num_neg].tolist()
    y_red = np.ones(int(num_pos+num_neg))*(-1.)
    y_red[0:num_pos]=1.
#y_red=y_red.astype(int)
    return x_red,y_red

def svmTraining(x,y,C,kernel_degree):
    feats_train  = StringCharFeatures(x,DNA)
    labels = BinaryLabels(y);
    start = timeit.default_timer()
    print "compute kernel matrix"
    kernel = WeightedDegreePositionStringKernel(feats_train, feats_train, kernel_degree) 
    stop = timeit.default_timer()
    time_kernel = stop-start
    start = timeit.default_timer()
    svm = LibSVM(C, kernel, labels)
    print "train support vector machine"
    svm.train()
    stop=timeit.default_timer()
    time_svm = stop-start
    return [svm,time_kernel,time_svm]

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
    labels = BinaryLabels(np.array(y));
    print "compute kernel matrix"
    kernel = WeightedDegreePositionStringKernel(feats_train, feats_train, kernel_degree)
    C=1
    svm = LibSVM(C, kernel,labels)
    print "train support vector machine"
    svm.train()
    pdb.set_trace()
    tally=len(x[0])
    print "compute poim"
    ma = poim.compute_poims(svm,kernel,poim_degree,tally)
    fobj = open(savepath,'wb')
    pickle.dump(ma,fobj)
    fobj.close()    
    return ma
