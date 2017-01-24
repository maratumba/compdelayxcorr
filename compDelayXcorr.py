#!/usr/bin/python
# cross correlate seismogram pairs and print out delays
# 
# awk 'NF<6 {split($1,a,".");split($2,b,".");print a[3],b[3],$5,$0}' xcorDelays4.txt > xcorDelays4Formatted.txt
import argparse
from obspy.taup.tau import TauPyModel

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Calculates differential arrival times')
    parser.add_argument('-e','--evinfo',type=str,help='file with event information',required=True)
    parser.add_argument('-o','--outfile',type=str,help='write output to a file')
    parser.add_argument('-v','--verbose',action='store_true',help='verbose output')
    parser.add_argument('FILES',nargs='+')
    args=parser.parse_args()

import glob
import obspy
import obspy.signal.cross_correlation
from obspy.core import read
import scipy.signal
import numpy as np
from obspy.taup.taup import getTravelTimes
from obspy.core.utcdatetime import UTCDateTime
import os
import sys
import pickle
  
def haversine(lat1,lon1,lat2,lon2):
    """
    haversine(lat1,lon1,lat2,lon2):
    Return distance in kilometers
    """
    
    
    import math
    R = 6372.8
    # In kilometers
    dLat = math.radians(lat2 - lat1)
    dLon = math.radians(lon2 - lon1)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)
    
    a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.sin(dLon / 2) * math.sin(dLon / 2) * math.cos(lat1) * math.cos(lat2)
    c = 2 * math.asin(math.sqrt(a))
    return R * c

def haversineD(lat1,lon1,lat2,lon2):
    """
    haversine(lat1,lon1,lat2,lon2):
    Return distance in degrees
    """
    
    
    import math
    R = 6372.8
    # In kilometers
    dLat = math.radians(lat2 - lat1)
    dLon = math.radians(lon2 - lon1)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)
    
    a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.sin(dLon / 2) * math.sin(dLon / 2) * math.cos(lat1) * math.cos(lat2)
    c = 2 * math.asin(math.sqrt(a))
    return c * 180 / math.pi

def parseEventDb(eventDbFile):
  #parse event catalog into a dictionary 
  colnames=['year','julday','time','depth','evla','evlo','evmag','evmag2']
  dtypes=['i','i','S8','f','f','f','f','f']
  dtype=zip(colnames,dtypes)
  evdb=np.genfromtxt(eventDbFile, dtype=dtype,)
  evdict={}
  
  for ev in evdb:
    dirname="_".join([str(ev['year']),"%03d" % ev['julday'],ev['time'],"M"+"%2.1f"%ev['evmag']])
    evdict[dirname]=ev
  
  
  return evdict

if __name__=='__main__':

    #sampling rate for cross correlation
    #crossDelta=0.001
    crossDelta=0.005
    
    #before and after wrt arrival time in seconds
    before=10.0
    after=20.0
    
    
    #skip counters, filename as keys
    #no travel time
    n_badtt={}
    #no distance
    n_badgcarc={}
    #no data in window
    n_badwin={}
    #no sync
    n_nosync={}
    #no BAZ
    n_nobaz={}
    
    
    #limit for difference in cross correlation waveforms start times
    sync_tolerance=0.01 #secs
    
    toolbar_width=40
    #cross correlation shift length
    #1/2 * window length
    #shiftlen=7500
    shiftlen=int(((before+after)/crossDelta)/4)
    
    stcut=obspy.core.stream.Stream()   
    #delays list
    delays=[]
    
    delaydict={}
    
    #parse event catalog
    evdict=parseEventDb(args.evinfo)
    dirname=os.getcwd().split('/')[-1]
    try:
        ev=evdict[dirname]
    except:
        print "ERROR: event "+dirname+" is not in the database, quitting"
        sys.exit(1)
    
    evyear=int(ev['year'])
    evjulday=int(ev['julday'])
    evhour=int(ev['time'][0:2])
    evmin=int(ev['time'][2:4])
    evsec=int(ev['time'][4:6])
    evmsec=int(ev['time'][7:8])*1000
    evtime=UTCDateTime(year=evyear,julday=evjulday,hour=evhour,minute=evmin,second=evsec,microsecond=evmsec)
    
    
    #initialize progress bar:
    nfiles=len(args.FILES)**2
    iprocprev=0
    sys.stdout.write("Progress: [%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1))
    ifile=0
    
    #create travel time calculator object
    taupmodel=TauPyModel(model='ak135')
    
    for bhz1File in args.FILES:
      if args.verbose:
          print "reading",bhz1File
      
      st1=read(bhz1File)
      tr1=st1[0]
      
      delta=tr1.stats.delta
      #remove mean:
      tr1.data=tr1.data-tr1.stats.sac.depmen
      
      #remove data before B sac header value
      #tr1.data=tr1.data[-tr1.stats.sac.b/tr1.stats.delta:]
      try:
        #calculate GCARC (for some reason obspy does not read gcarc header)
        gcarc1=haversineD(tr1.stats.sac.evla, tr1.stats.sac.evlo, tr1.stats.sac.stla, tr1.stats.sac.stlo)
      except:
        print "ERROR: Can\'t calculate event distance, quitting",tr1
        #unlikely that any other station will have a P arrival, so break
        sys.exit(1)
      
        
        print "delta=",gcarc1,"depth=",tr1.stats.sac.evdp
      #t1 = np.arange(0, tr1.stats.npts / tr1.stats.sampling_rate, tr1.stats.delta)
      
      
      try:
        #tt1=getTravelTimes(delta=gcarc1, depth=tr1.stats.sac.evdp, model='ak135',phase_list=['P'])[0]['time']
        tt1=taupmodel.get_travel_times(tr1.stats.sac.evdp, gcarc1, ['P'])[0].time
      except:
        print "WARNING: Can\'t calculate travel time, skipping (No P arrival for distance:",gcarc1,")",bhz1File
        continue
      

      tr1.stats.sac.t1=evtime+tt1
      tr1.trim(tr1.stats.sac.t1-before, tr1.stats.sac.t1+after)
      if len(tr1.data)==0:
          if args.verbose:
              print "WARNING: no data within window, skipping",bhz1File
          continue
          
      
      try:
          crossNpts=len(tr1.data)/(crossDelta/tr1.stats.delta)
          #sys.exit()
          
          crossData1=scipy.signal.resample(tr1.data,crossNpts)
      except:
          print "Problem resampling, skipping",bhz1File
          continue
          
          
      if tr1.stats.station not in delaydict:
          delaydict[tr1.stats.station]={}
      
      for bhz2File in args.FILES:
        ifile+=1
        if args.verbose:
            print "reading",bhz2File
        else:
            iproc=int(float(ifile)/nfiles*40.)
            if iproc>iprocprev:
                sys.stdout.write('#'*(iproc-iprocprev))
                sys.stdout.flush()
                iprocprev=iproc
                   
        st2=read(bhz2File)
        tr2=st2[0]  
        
        try:
          #calculate GCARC (for some reason obspy does not read gcarc header)
          gcarc2=haversineD(tr1.stats.sac.evla, tr1.stats.sac.evlo, tr2.stats.sac.stla, tr2.stats.sac.stlo)
        except:
          print "WARNING: Can\'t calculate event distance, skipping",tr2
          continue
          
          print "delta=",gcarc2,"depth=",tr2.stats.sac.evdp
        try:
          #tt2=getTravelTimes(delta=gcarc2, depth=tr2.stats.sac.evdp, model='ak135',phase_list=['P'])[0]['time']
          tt2=taupmodel.get_travel_times(tr2.stats.sac.evdp, gcarc2, ['P'])[0].time
        except:
          print "WARNING: Can\'t calculate travel time, skipping (No P arrival for distance:",gcarc2,")",bhz2File
          continue
        
        tr2.stats.sac.t1=evtime+tt2
        tr2.data=tr2.data-tr2.stats.sac.depmen
        #tr2.data=tr2.data[-tr2.stats.sac.b/tr2.stats.delta:]
        tr2.trim(tr1.stats.sac.t1-before, tr1.stats.sac.t1+after)
        
        if len(tr2.data)==0:
            if args.verbose:
                print "WARNING: no data within window, skipping",bhz2File
            continue
                #make sure start times are equal
        if abs(tr1.stats.starttime-tr2.stats.starttime)>sync_tolerance:
          #I can fix this
          print "WARNING: starttimes diff>"+sync_tolerance+" skipping",bhz1File,bhz2File,tr1.stats.starttime,tr1.stats.starttime-tr2.stats.starttime
          continue 
        
        try:
            crossData2=scipy.signal.resample(tr2.data,crossNpts)
        except:
            print "Problem resampling, skipping",bhz2File
            sys.exit()
            continue
        
        maxInd,maxval=obspy.signal.cross_correlation.xcorr(crossData1,crossData2,shiftlen)
        
        t2del=maxInd*crossDelta
        t1del=tr1.stats.sac.t1 - tr2.stats.sac.t1
        
        deldel=t2del-t1del
        #delays.append([tr1.__str__(),tr2.__str__(),tr1.stats.station,tr2.stats.station,t2del,t1del,deldel,maxval])
        #stcut.append(tr2)
        #pairstr=tr1.stats.station+'_'+tr2.stats.station
        
        if tr2.stats.station not in delaydict[tr1.stats.station]:
            delaydict[tr1.stats.station][tr2.stats.station]={}
            delaydict[tr1.stats.station][tr2.stats.station]['t1del']=t1del    
            delaydict[tr1.stats.station][tr2.stats.station]['t2del']=t2del
            delaydict[tr1.stats.station][tr2.stats.station]['deldel']=deldel
            delaydict[tr1.stats.station][tr2.stats.station]['maxval']=maxval
            delaydict[tr1.stats.station][tr2.stats.station]['sac1']=bhz1File
            delaydict[tr1.stats.station][tr2.stats.station]['sac2']=bhz2File
            delaydict[tr1.stats.station][tr2.stats.station]['starttime']=tr1.stats.starttime
            delaydict[tr1.stats.station][tr2.stats.station]['endtime']=tr1.stats.endtime
            try:
                delaydict[tr1.stats.station][tr2.stats.station]['baz']=tr1.stats.sac.baz
            except:
                if args.verbose:
                    print "WARNING: No BAZ in SAC file",bhz1File
                delaydict[tr1.stats.station][tr2.stats.station]['baz']=-12345
            
            delaydict[tr1.stats.station][tr2.stats.station]['gcarc1']=gcarc1
        else:
            print "WARNING: station pair",tr1.stats.station,tr2.stats.station,"already in record"
            
      
        if args.verbose:
            print bhz1File, bhz2File, t2del, t1del, deldel, maxval
        #print tr1.stats.station, tr2.stats.station, t2del, t1del, deldel
    sys.stdout.write(']\n')
    sys.stdout.flush()
    
    pickle.dump(delaydict,open(dirname+'.pickle','wb'))
    
    #convert dict to numpy structured array
    with open(dirname+'.txt','w') as f:
        delaylist=[]
        for st1 in delaydict:
            for st2 in delaydict[st1]:
                line=[]
                line.append(st1)
                line.append(st2)
                line.append(delaydict[st1][st2]['t1del'])
                line.append(delaydict[st1][st2]['t2del'])
                line.append(delaydict[st1][st2]['deldel'])
                line.append(delaydict[st1][st2]['maxval'])
                line.append(delaydict[st1][st2]['sac1'])
                line.append(delaydict[st1][st2]['sac2'])
                line.append(delaydict[st1][st2]['starttime'])
                line.append(delaydict[st1][st2]['endtime'])
                line.append(delaydict[st1][st2]['baz'])
                #print st1,st2,delaydict[st1][st2][valname]
                    #break
                f.write(" ".join(map(str,line))+"\n")
                delaylist.append(line)
                
        
    
    
    
    
    
    
    
