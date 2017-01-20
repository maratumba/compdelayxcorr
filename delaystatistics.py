#!/usr/bin/env python

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_delays(delaydict):
    
    deldels=delaydict['deldel']
    mean=np.mean(deldels)
    std=np.std(deldels)
    
    xmax=len(deldels)+1
    xmin=0
    ymax=1.0
    ymin=-1.0
    
    #plot uncertainty bar
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    
    ax1.set_ylim(-1,1)
    ax1.set_xlim(0,xmax)
    
    ax1.add_patch(
        patches.Rectangle(
            (xmin, mean-std),   # (x,y)
            xmax,          # width
            2*std,          # height
            alpha=0.3,
            color='blue',
        )
    )
    
    #plot the mean line
    ax1.plot([xmin,xmax],[mean,mean],'k-',alpha=0.5)
    
    ax1.plot(range(1,len(deldels)+1),deldels,'o')

    
    plt.show()

def plot_delays_ensemble(delaydict):
    pass


def shed_outliers(delaydict,col='deldel',std_lim=1,n_iter=2):
    '''
    Returns an array without the outliers.
    std_lim -> threshold as std factor
    n_iter -> number of iterations
    col -> name of the colon in the input numpy array 
    '''
    
    
    d=delaydict.copy()
    
    for step in range(n_iter):
        mean=np.mean(d[col])
        std=np.std(d[col])
        d2=np.array([],dtype=d.dtype)
        for r in d:
            if r[col]>mean-std*std_lim and r[col]<mean+std*std_lim:
                d2=np.append(d2,np.array(r,dtype=d.dtype))
        d=d2
    
    return d
            
    
    
    
    

if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description='Calculates differential arrival times')
    parser.add_argument('-o','--outfile',type=str,help='write output to a file')
    parser.add_argument('-v','--verbose',action='store_true',help='verbose output')
    parser.add_argument('FILES',nargs='+')
    args=parser.parse_args()
    
    
    delaydict={}
    
    colnames=['t1del','t2del','deldel','maxval']
    dtypes=['f4','f4','f4','f4']
    dtype=zip(colnames,dtypes)
    
    #parse and gather txt files
    for txt_file in args.FILES:
        with open(txt_file,'r') as f:
            for line in f.readlines():
                line=line.split()
                if line[0]+"_"+line[1] not in delaydict:
                    delaydict[line[0]+"_"+line[1]]=np.array([],dtype=dtype)
                #if line[1] not in delaydict[line[0]]:
                #    delaydict[line[0]][line[1]]=np.array([],dtype=dtype)
                
                deldel=np.array(tuple(line[2:6]),dtype=dtype)
                delaydict[line[0]+"_"+line[1]]=np.append(delaydict[line[0]+"_"+line[1]],deldel)
           

    
             
    #for st1 in delaydict:
    #    for st2 in delaydict[st1]:
    #        print delaydict[st1][st2]
                
