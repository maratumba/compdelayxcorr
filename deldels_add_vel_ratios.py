#!/usr/bin/env python

import sys
import argparse
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import re
from numpy.lib.recfunctions import append_fields

def plot_v12s(delaydict,key,outfile=None):

    deldels=delaydict['v12']
    mean=np.mean(deldels)
    std=np.std(deldels)

    xmax=len(deldels)+1
    xmin=0
    ymax=1.5
    ymin=0.5

    #plot uncertainty bar
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    ax1.set_ylim(ymin,ymax)
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

def shed_outliers_list(vlist,std_lim=1,n_iter=2):
    '''
    Returns an array without the outliers.
    std_lim -> threshold as std factor
    n_iter -> number of iterations
    col -> name of the colon in the input numpy array
    '''
    d=np.array(vlist).copy()

    #print np.mean(d['v'][col])

    for step in range(n_iter):
        #print step
        mean=np.mean(d)
        std=np.std(d)
        d2=[]
        for r in d:
            #print(r)
            if r>mean-std*std_lim and r<mean+std*std_lim and std!=0:
                #print r,'included'
                d2.append(r)
            else:
                #print r,'excluded'
                pass
       #print d
        #print d2
        mean=np.mean(d2)
        std=np.std(d2)

        d=d2

    return d,mean,std

def get_content(re1,dict1):

    cont={}
    for key in dict1:
        if re.match(re1,key):
            cont[key]=dict1[key]
    return cont

def parseStationDb(stationDbFile):
    colnames=['stname','stlat','stlon','stel']
    dtypes=['S8','f','f','f']
    dtype=zip(colnames,dtypes)
    stdb=np.genfromtxt(stationDbFile, dtype=dtype,)
    stdict={}

    for st in stdb:
        stdict[st['stname']]=st

    return stdict

def calc_vel_ratio(delaydict, h=35, v2=6.3):
    """
    v1/v2=h/(h+deldel*v2*cos(theta))
    :param delaydict:
    :return:
    """
    # h=35 #km
    # v2=6.3 #km/s
    colnames=['t1del','t2del','deldel','maxval','sac1','sac2','starttime','endtime','ar1','ar2','baz','incident_angle_rad_1','incident_angle_rad_2','v12']
    dtypes=['f4','f4','f4','f4','S40','S40','S40','S40','S40','S40','f4','f4','f4','f4']
    dtype=zip(colnames,dtypes)
    delaydict2={}
    for pair in delaydict:
        l=len(delaydict[pair]['v'])
        a=3
        v12_col=np.zeros(l, dtype=[('v12', 'f4')])
        #if 'v12' not in delaydict[pair]['v'].dtype.names:
        # delaydict[pair]['v']=append_fields(delaydict[pair]['v'], 'v12', v12_col, usemask=False)
        delaydict2[pair]={}
        delaydict2[pair]['v']=np.array([],dtype=dtype)
        for v in delaydict[pair]['v']:
            v12=h/(h+v['deldel']*v2*np.cos(v['incident_angle_rad_1']))
            newline=np.array(tuple(list(v)+[v12]),dtype=dtype)
            delaydict2[pair]['v']=np.append(delaydict2[pair]['v'],newline)
            #import pdb;pdb.set_trace()
            #v['v12']=v12

    return delaydict2

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


def shed_bad_corr(delaydict,col='deldel',corr_lim=0.6):
    '''
    Returns an array without the outliers.
    std_lim -> threshold as std factor
    n_iter -> number of iterations
    col -> name of the colon in the input numpy array
    '''

    d=delaydict.copy()

    #print np.mean(d['v'][col])
    #import pdb; pdb.set_trace()

    d2={'v':np.array([],dtype=d['v'].dtype)}

    for r in d['v']:
        if r['maxval']>corr_lim:
            d2['v']=np.append(d2['v'],np.array(r,dtype=d['v'].dtype))


    d2['mean']=np.mean(d2['v'][col])
    d2['std']=np.std(d2['v'][col])

    d['v']=d2['v']
    d['mean']=d2['mean']
    d['std']=d2['std']

    return d


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


def plot_profiles_ref(delaydict,stdict,outfile=None,ymin=0.5,ymax=1.5):

    mpl.style.use('ggplot')

    fig=plt.figure()
    i=1
    for c in range(ord('A'),ord('F')+1):
        keys=get_content(r'D'+chr(c)+r'01_D'+chr(c)+r'.*',delaydict).keys()

        means=[]
        for stpair in keys:
            print stpair
            #import pdb;pdb.set_trace()
            means.append(shed_outliers(shed_bad_corr(delaydict[stpair],corr_lim=0.7,col='v12'),col='v12',std_lim=1,n_iter=1)['mean'])


        #means=[shed_outliers(shed_bad_corr(delaydict[x],corr_lim=0.7,col='v12'),col='v12',std_lim=4,n_iter=1)['mean'] for x in keys]
        stds=[shed_outliers(shed_bad_corr(delaydict[x],corr_lim=0.7,col='v12'),col='v12',std_lim=1,n_iter=1)['std'] for x in keys]

        lats=[stdict[x[5:9]]['stla'] for x in keys]

        means=[x for (y,x) in sorted(zip(keys,means))]
        stds=[x for (y,x) in sorted(zip(keys,stds))]
        lats=[x for (y,x) in sorted(zip(keys,lats))]


        ax=plt.subplot(6,1,i)
        baseline,=plt.plot(lats,means)#,label='D'+chr(c))
        ticks=plt.plot(lats,means,'o')
        lc=baseline.get_color()
        plt.fill_between(lats, np.array(means)-np.array(stds), np.array(means)+np.array(stds),alpha=0.2,facecolor=lc)

        ax.set_ylim(ymin, ymax)

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

        plt.yticks([ymin,1-(1-ymin)/2,1,1+(ymax-1)/2,ymax],["%0.2f" % ymin,'','1.0','',"%0.2f" % ymax])

        plt.legend(frameon=False)
        i+=1

    if outfile:
        fig.savefig(outfile)
    else:
        plt.show()



if __name__=='__main__':

    parser=argparse.ArgumentParser(description='Calculates differential arrival times statistics')
    parser.add_argument('-o','--outfile',type=str,help='write output to a file')
    parser.add_argument('-v','--verbose',action='store_true',help='verbose output')
    parser.add_argument('FILES',nargs='+')
    args=parser.parse_args()

    delaydict_file = args.FILES[0]
    delaydict_file_pre = '.'.join(delaydict_file.split('.')[:-1])
    delaydict = pickle.load(open(delaydict_file,'rb'))
    # calc_vel_ratio(delaydict)
    stdict = parse_station_info('station_info.txt')
    #plot_profiles_ref(delaydict,stdict)

    delaydict_35__6_3 = calc_vel_ratio(delaydict, h=35, v2=6.3)
    delaydict_10__6_3 = calc_vel_ratio(delaydict, h=10, v2=6.3)

    pickle.dump(delaydict_35__6_3, open(delaydict_file_pre+'_v12_35__6_3.pickle', 'wb'))
    pickle.dump(delaydict_10__6_3, open(delaydict_file_pre+'_v12_10__6_3.pickle', 'wb'))

    print_grid_ref(delaydict_35__6_3, stdict, outfile=delaydict_file_pre + '_v12_35__6_3.txt')
    print_grid_ref(delaydict_10__6_3, stdict, outfile=delaydict_file_pre + '_v12_10__6_3.txt')

    plot_profiles_ref(delaydict_35__6_3, stdict, outfile=delaydict_file_pre + '_v12_35__6_3.png')
    plot_profiles_ref(delaydict_10__6_3, stdict, outfile=delaydict_file_pre + '_v12_10__6_3.png')






