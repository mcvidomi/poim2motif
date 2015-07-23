'''
Created on 16.06.2015

@author: marinavidovic
'''

import numpy as np
import math
import matplotlib
import pdb
import matplotlib
matplotlib.use('Agg')
import pylab as py
def test():

    plt.plot([1,2,3])
    plt.savefig('testBB.png')


def figurepoimsimple(poim, savefile, show):
    matplotlib.use('Agg')
    #py.figure(figsize=(14, 12))
    py.ioff()
    motivelen = int(np.log(len(poim)) / np.log(4))
    ylabel = []
    for i in range(int(math.pow(4, motivelen))):
        label = []
        index = i
        for j in range(motivelen):
            label.append(index % 4)
            index = int(index / 4)
        label.reverse()
        ylabel.append(veclisttodna(label))
        
    py.pcolor(poim)
    
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(40)  
    py.axis([0, len(poim[0]), 0, len(poim)])
    diff = int((len(poim[0]) / 5)) - 1
    x_places = py.arange(4.5, len(poim[0]), diff)
    py.xticks(x_places, np.arange(5, len(poim[0]) + 1, diff),fontsize=40) 
    
    py.xlabel("Position", fontsize=46)
    py.ylabel("Motif", fontsize=46)
    py.yticks(np.arange(math.pow(4, motivelen)) + 0.5, (ylabel), fontsize=40)   
    if savefile != "":
#        F = py.gcf()
#        # Now check everything with the defaults:
#        DefaultSize = F.get_size_inches()
#        # Now make the image twice as big, while keeping the fonts and all the
#        # same size
#        F.set_size_inches((DefaultSize[0]*2, DefaultSize[1]*2))
        py.savefig(savefile) 
    if show:
        py.show()
    
def figurepoimsimple_small(poim, l, start, savefile, show):
    R = poim

    py.figure(figsize=(14, 12))
    motivelen = int(np.log(len(poim)) / np.log(4))
    ylabel = []
    for i in range(int(math.pow(4, motivelen))):
        label = []
        index = i
        for j in range(motivelen):
            label.append(index % 4)
            index = int(index / 4)
        label.reverse()
        ylabel.append(veclisttodna(label))
    py.pcolor(R[:, start:start + l])
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
     t.set_fontsize(40)   
    diff = int((l / 5)) - 1
    x_places = py.arange(0.5, l, diff)
    xa = np.arange(start, start + l, diff)
    diff = int((l / 4)) 
    x_places = py.arange(0.5, l , diff)
    xa = np.arange(start + 1, start + 1 + l, diff)
    py.xlabel("Position", fontsize=46)
    py.ylabel("Motif", fontsize=46)
    py.yticks(np.arange(math.pow(4, motivelen)) + 0.5, (ylabel),fontsize=40)   
    py.xticks(x_places, (xa.tolist()),fontsize=40)
    if savefile != "":
        py.savefig(savefile)  
    print "the poim should show up here"
    if show:
        py.show()


def figuremaxPOIM_small(poim, savepath, l , start, show):  
    R = poim

    py.figure(figsize=(14, 12))
    py.pcolor(R[:, start:start + l])
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
     t.set_fontsize(40)
    diff = int((l / 4)) 
    x_places = py.arange(0.5, l , diff)
    xa = np.arange(start + 1, start + 1 + l, diff)
    py.xlabel("Position", fontsize=46)
    py.ylabel("Motif length", fontsize=46)
    py.xticks(x_places, (xa.tolist()),fontsize=40)
    py.yticks(np.arange(0, len(poim)) + 0.5, np.arange(1, len(poim) + 1),fontsize=40)    
    py.savefig(savepath)
    if show:
        py.show()

  
def figuremaxPOIM(R, savepath, show):  

    py.figure(figsize=(14, 12))
    py.pcolor(R)
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(40)
    py.axis([0, len(R[0]), 0, len(R)])
    diff = int((len(R[0]) / 5)) - 1
    x_places = py.arange(4.5, len(R[0]), diff)
    py.xticks(x_places, np.arange(5, len(R[0]) + 1, diff),fontsize=40) 
    py.xlabel("Position", fontsize=46)
    py.ylabel("Motif length", fontsize=46)
    py.yticks(np.arange(0, len(R)) + 0.5, np.arange(1, len(R) + 1),fontsize=40)    
    py.savefig(savepath)
    if show:
        py.show()

        
    
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