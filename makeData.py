import random
import cPickle as pickle
from numpy import concatenate, ones, array, shape, size, zeros, exp, arange
from numpy import concatenate,ones,array,shape,size,zeros,exp
import numpy as np
import copy
import math
import pdb
dna = ['A', 'C', 'G', 'T']

def simulate_sequence(length):
    sequence = ''
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


def gensequences2(tally,positives,sequenceno,prob,motif,mu):
    sequences = []
    ml = len(motif)
    for i in range(sequenceno):
        aa = simulate_sequence(tally)
        if i < positives:
            
            mut=mutate_motif(motif,prob)
            aa = aa.replace(aa[mu:mu + ml], mut)
        sequences.append(aa)
    return sequences

def gensequences(tally,positives,sequenceno,prob,motif,mu):
    sequences = []
    y=np.ones(sequenceno)*(-1)
    ml = len(motif)
    for i in range(sequenceno):
        aa = simulate_sequence(tally)
        if i < positives:
            y[i]=1
            mut=mutate_motif(motif,prob)
            aa = aa.replace(aa[mu:mu + ml], mut)
        sequences.append(aa)
    return sequences,y

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
    for i in range(len(counter[0])):
        for j in range(4):
            if counter[j,i] == 1.0:
                print "nucleotid " , dna[j]," position", str(i)  
    

    return counter

def extractRealData(datapath,savepath,lines):
#path = "/home/mvidovic/POIMall/data/real/"
#filename = "human_acceptor_splice_data.txt"
    data =file(datapath).readlines()[:lines]
    labels = []
    x = []
    cn=0
    for i in range(len(data)):
        labels.append(int(data[i][0:2]))
        x.append(data[i][3:-1])
        pdb.set_trace()
        
        
        if int(data[i][0:2])==1:
            cn=cn+1
    print "numper of positive labels: " , cn
    if savepath !="":
        fobj = open(savepath,'wb')
        pickle.dump([x,labels],fobj) 
        fobj.close()
    return x,labels




#tally = 30       #length of training sequences
#sequenceno = 50 #number of training sequences
#positives = 12 #positives training sequences
#motiv = "CCTATA"
#mu = 10
#poim_degree = 3
#idxprob=np.arange(0.0,1.01,1.1)
#dna = ['A', 'C', 'G', 'T']

def compute_data(seq_length,seq_no,pos_no,motifs,mu):
    a=1

#fm_train_dna = gensequences(tally,positives,sequenceno,0,motiv,mu)
