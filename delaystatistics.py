#!/usr/bin/env python

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import re
import pickle
from obspy.core import read
from obspy.core.utcdatetime import UTCDateTime
import glob

def plot_shift_wf(delaydict,basedir='',ax1=None):

    #delaydict=delaydict[1]

    shift=delaydict['v']['t2del']
    cc=delaydict['v']['maxval']

    trim_left=5.8 #12.0
    trim_right=19.2 #15.6

    sac1=glob.glob(basedir+'/*/'+delaydict['v']['sac1'])[0]
    sac2=glob.glob(basedir+'/*/'+delaydict['v']['sac2'])[0]

    tr1=read(sac1)[0]
    tr2=read(sac2)[0]
    tr2_o=tr2.copy()

    tr1_t1=UTCDateTime(delaydict['v']['starttime'])+trim_left
    tr1_t2=UTCDateTime(delaydict['v']['endtime'])-trim_right

    tr2_t1=tr1_t1-shift
    tr2_t2=tr1_t2-shift

    tr1.data=tr1.data-np.mean(tr1.data[0:100])
    tr2.data=tr2.data-np.mean(tr2.data[0:100])

    tr1.trim(tr1_t1,tr1_t2)
    tr2.trim(tr2_t1,tr2_t2)
    tr2_o.trim(tr1_t1,tr1_t2)

    t=tr1.times()


    if not ax1:
        with plt.style.context('ggplot'):
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)

    with plt.style.context('ggplot'):

        ax1.plot(t,tr1.normalize().data,label='tr1')
        ax1.plot(t,tr2.normalize().data,label='tr2')
        #ax1.plot(t,tr2_o.normalize().data,label='tr2_o')

        ax1.text(0.01, 0.99,'shift='+str(shift)+'\ncc='+str(cc),
                 verticalalignment='top', horizontalalignment='left',
                 transform=ax1.transAxes,
                 color='black', fontsize=15)


        ax1.tick_params(
                        which='both',      # both major and minor ticks are affected
                        bottom='off',      # ticks along the bottom edge are off
                        top='off',         # ticks along the top edge are off
                        left='off',
                        right='off')
        plt.legend(frameon=False)
        plt.show()


def get_content(re1,dict1):

    cont={}
    for key in dict1:
        if re.match(re1,key):
            cont[key]=dict1[key]
    return cont

def plot_delays(delaydict,key,outfile=None):

    deldels=delaydict['v']['deldel']
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
    ax1.set_title(key.split('_')[0]+'-'+key.split('_')[1])
    ax1.text(0.5, 0.95,'mean='+str(mean)+'\nstd='+str(std),
             verticalalignment='top', horizontalalignment='center',
             transform=ax1.transAxes,
             color='green', fontsize=15)
    if outfile:
        fig1.savefig(outfile)
    else:
        plt.show()

def plot_delays_ensemble(delaydict):
    pass

def shed_bad_corr(delaydict,corr_lim=0.6):
    '''
    Returns an array without the outliers.
    std_lim -> threshold as std factor
    n_iter -> number of iterations
    col -> name of the colon in the input numpy array
    '''

    d=delaydict.copy()

    #print np.mean(d['v'][col])


    d2={'v':np.array([],dtype=d['v'].dtype)}

    for r in d['v']:
        if r['maxval']>corr_lim:
            d2['v']=np.append(d2['v'],np.array(r,dtype=d['v'].dtype))

    d2['mean']=np.mean(d2['v']['deldel'])
    d2['std']=np.std(d2['v']['deldel'])

    d['v']=d2['v']
    d['mean']=d2['mean']
    d['std']=d2['std']

    return d

def shed_outliers(delaydict,col='deldel',std_lim=1,n_iter=2):
    '''
    Returns an array without the outliers.
    std_lim -> threshold as std factor
    n_iter -> number of iterations
    col -> name of the colon in the input numpy array
    '''


    d=delaydict.copy()

    #print np.mean(d['v'][col])

    for step in range(n_iter):
        #print step
        mean=np.mean(d['v'][col])
        std=np.std(d['v'][col])
        d2={'v':np.array([],dtype=d['v'].dtype)}
        for r in d['v']:
            if r[col]>mean-std*std_lim and r[col]<mean+std*std_lim and std!=0:
                d2['v']=np.append(d2['v'],np.array(r,dtype=d['v'].dtype))
       #print d
        #print d2
        d2['mean']=np.mean(d2['v'][col])
        d2['std']=np.std(d2['v'][col])

        d['v']=d2['v']
        d['mean']=d2['mean']
        d['std']=d2['std']

    return d

def parse_station_info(stinfo_file):
    stdict={}
    colnames=['stla','stlo','stel']
    dtypes=['f4','f4','f4']
    dtype=zip(colnames,dtypes)

    with open(stinfo_file) as f:
        for line in f.readlines():
            line=line.split()
            stinfo=np.array(tuple(line[1:4]),dtype=dtype)
            stdict[line[0]]=stinfo

    return stdict

def print_grid_ref(delaydict,stdict,refst='DA01',outfile=None,do_shed_outliers=True):


    keys=get_content(refst+r'_.*',delaydict).keys()

    if outfile:
        of=open(outfile,'w')

    for key in delaydict.keys():
        if do_shed_outliers:
            delaydict[key]=shed_outliers(delaydict[key])

        mean=delaydict[key]['mean']
        stname=key.split('_')[1]
        outstr=' '.join([stname,str(stdict[stname]['stlo']),str(stdict[stname]['stla']),str(mean)])
        if not outfile:
            print outstr
        else:
            of.write(outstr+'\n')

    if outfile:
        of.close()


def plot_profiles_ref(delaydict,stdict):

    mpl.style.use('ggplot')

    fig=plt.figure()
    i=1
    for c in range(ord('A'),ord('F')+1):
        keys=get_content(r'D'+chr(c)+r'01_D'+chr(c)+r'.*',delaydict).keys()

        means=[shed_outliers(shed_bad_corr(delaydict[x],corr_lim=0.7),std_lim=4,n_iter=1)['mean'] for x in keys]
        stds=[shed_outliers(shed_bad_corr(delaydict[x],corr_lim=0.7),std_lim=4,n_iter=1)['std'] for x in keys]

        lats=[stdict[x[5:9]]['stla'] for x in keys]

        means=[x for (y,x) in sorted(zip(keys,means))]
        stds=[x for (y,x) in sorted(zip(keys,stds))]
        lats=[x for (y,x) in sorted(zip(keys,lats))]


        ax=plt.subplot(6,1,i)
        baseline,=plt.plot(lats,means)#,label='D'+chr(c))
        ticks=plt.plot(lats,means,'o')
        lc=baseline.get_color()
        plt.fill_between(lats, np.array(means)-np.array(stds), np.array(means)+np.array(stds),alpha=0.2,facecolor=lc)


        ax.set_ylim(-0.4,0.4)

        plt.tick_params(
                    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    left='off',
                    right='off')

        plt.suptitle('Deldel wrt station 01 in each line')

        ax.text(0.99, 0.95,'D'+chr(c),
                verticalalignment='top', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=15)

        #ax.set_yticks([-0.4,-0.2,0,0.2,0.4],['-0.4','','0.0','','0.4'])

        if i!=6:
            ax.axes.get_xaxis().set_ticklabels([])

        if i==6:
            plt.xlabel('Latitude (deg)')
            plt.ylabel('Mean deldel (s)')

        plt.yticks([-0.4,-0.2,0,0.2,0.4],['-0.4','','0.0','','0.4'])

        plt.legend(frameon=False)
        i+=1

    plt.show()

if __name__=='__main__':

    parser=argparse.ArgumentParser(description='Calculates differential arrival times statistics')
    parser.add_argument('-o','--outfile',type=str,help='write output to a file')
    parser.add_argument('-v','--verbose',action='store_true',help='verbose output')
    parser.add_argument('FILES',nargs='+')
    args=parser.parse_args()


    delaydict={}

    colnames=['t1del','t2del','deldel','maxval','sac1','sac2','starttime','endtime','ar1','ar2','baz','incident_angle_rad_1','incident_angle_rad_2']
    dtypes=['f4','f4','f4','f4','S40','S40','S40','S40','S40','S40','f4','f4','f4']
    dtype=zip(colnames,dtypes)

    #parse and gather txt files
    for txt_file in args.FILES:
        print('Reading {}'.format(txt_file))
        with open(txt_file,'r') as f:
            for line in f.readlines():
                line=line.split()
                try:
                    if line[0]+"_"+line[1] not in delaydict:
                        delaydict[line[0]+"_"+line[1]]={'v':np.array([],dtype=dtype)}
                    #if line[1] not in delaydict[line[0]]:
                    #    delaydict[line[0]][line[1]]=np.array([],dtype=dtype)

                    # horrendous parsing:

                    deldel=np.array(tuple(line[2:]),dtype=dtype)
                    delaydict[line[0]+"_"+line[1]]['v']=np.append(delaydict[line[0]+"_"+line[1]]['v'],deldel)
                except:
                    print('Problem reading line {} from file {}'.format(line,txt_file))          


    for pair in delaydict:
        delaydict[pair]['mean']=np.mean(delaydict[pair]['v']['deldel'])
        delaydict[pair]['std']=np.std(delaydict[pair]['v']['deldel'])

    pickle.dump(delaydict,open('delays_collected.pickle','wb'))
    #stdict=parse_station_info('station_info.txt')
    #for st1 in delaydict:
    #    for st2 in delaydict[st1]:
    #        print delaydict[st1][st2]
