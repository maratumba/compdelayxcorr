#!/bin/bash

reg=-R29.0/31.0/40.0/41.5
proj=-JM6i

stLonLat=$1
gmtDataBase="/Users/gulyamani/research/gmtData"
xyz_maxFile=$1
ofile=$1.ps

echo wtf
#plot the points at max amp locations of the 3D stack between 20-51 km.

gmtset BASEMAP_TYPE fancy ANOT_FONT_SIZE 12
gmtset PLOT_DEGREE_FORMAT ddd.xx
gmtset PAPER_MEDIA letter+

psbasemap -P -X1i $reg $proj  -Ba0.50f1.0SWne -K > $ofile

makecpt -Cgray -Z -T-18000/4000/10 > grey.cpt
makecpt -Cpolar -T-0.5/0.5/0.01 -D -I -V > rainb.cpt
if [ ! -e ${gmtDataBase}/turkiye.grad ]; then grdgradient $gmtDataBase/turkiye.grd -G${gmtDataBase}/turkiye.grad -A -Nt0.4;fi

grdimage -R ${gmtDataBase}/turkiye.grd -JM -I${gmtDataBase}/turkiye.grad -P -K -O -Cgrey.cpt >> $ofile

faultfile=./faulttur.map

#block mean the maxima and fit a surface with high smoothness
#./xyzA2latlonzA.py $xyz_maxFile |\
#awk '$4>0.15 && $3>20 {print $2,$1,$3}' |\
awk '{print $2,$3,$4}' $xyz_maxFile |\
blockmean -fg -V $reg -I5k  |\
surface -fg -V $reg -I1 -T0.25 -Gdeldel.grd
#grdgradient soCalMoho.grd -G -A -Nt0.4
#grdimage $reg $proj -V -O -K -I${gmtDataBase}/turkiye.grad -Crainbo.cpt deldel.grd >> $ofile

#./xyzA2latlonzA.py $xyz_maxFile |\
#awk '$4>0.15 && $3>20 {print $2,$1,$3}' $xyz_maxFile |\
#awk '{print $2,$1,$3}' $xyz_maxFile |\
#blockmean -V $reg -I50k  > rf_maxZ_soCal.bmean.outF.txt

#coast lines, rivers and lakes
pscoast -V -O -K $reg $proj -W0.05 -Df -P -A100 -P -S255/255/255 >> $ofile
#pscoast -V -O -K $reg $proj -W0.05 -Df -P -A100 -P >> $ofile

#faults:
psxy -V -M -O -K $reg $proj -W3/100/100/100 $gmtDataBase/mta_emme_fault_map.gmt >> $ofile

#./xyzA2latlonzA.py $xyz_maxFile |\
#awk '$4>0.15 && $3>20 {print $2,$1,$3}' |\
awk '{print $2,$3,$4}' $xyz_maxFile |\
psxy -V -O -K $reg $proj -Crainb.cpt -Ss.3 -W0.1p >> $ofile

#scvm stations
#awk '{print $2,$3}' $xyz_maxFile |\
#psxy -V -O -K $reg $proj -G0/0/200 -Si.15 >> $ofile

#awk '{print $2,$3,"6 0 0 CM",$1}' $xyz_maxFile |\
#pstext -V -O -K -Y.1 $reg $proj >> $ofile

#legend
psscale -D16.6c/3.9c/8.3c/0.3c -E -B0.1/1:"mean deldel (s)": -O -K -Crainb.cpt -V >> $ofile

#finish him
echo 0 0 | psxy $reg $proj -Sa.1 -X.5i -O >> $ofile

echo writing $ofile

convert -density 300x300 $ofile $ofile.png
echo writing $ofile.png
