#!/bin/bash

#if ( "$#argv" < 7) 
#then
#  echo "Usage :  $0 -L[listFile] -D[RF dir] -P[vp] -Z[zmin/zmax/zinc]
#  -R[rmin/rmax/rinc] -W [w_2p1s/w_1p2s] -N[nth]"
#  echo "example: $0 -LRF.lst -DRF_dir -P6.3 -Z20/60/0.5 -R1.5/2.0/0.01
#  -W0.5/0.5 -N1"
#  exit
#fi

# Notice: 
#   1) zero time must be the zero-offset (P peak in P receiver function)
#   2) the file name should be less than 30 characters!

function usage(){
    printf "Usage: %s: [-L listFile] [-D RF_dir] [-P vp] [-Z zmin/zmax/zinc] [-R rmin/rmax/rinc] [-W w_2p1s/w_1p2s] [-N nth]\n" $0
    echo "Default: $0 [-LlistFile] -D. -P6.3 -Z20/60/0.5 -R1.5/2.0/0.01 -W0.5/0.5 -N1"
    exit 1
}
# Default parameters

staLIST=
RFdir=.
vp=6.3
zarg=20/60/0.5
rarg=1.5/2.0/0.01
w=0.5/0.5
nth=1

# parse options
while getopts L:D:P:Z:R:W:N: name
do
    case $name in
    L)  staLIST="$OPTARG";;
    D)  RFdir="$OPTARG";;
    P)  vp="$OPTARG";;
    Z)  zarg="$OPTARG";;
    R)  rarg="$OPTARG";;
    W)  w="$OPTARG";;
    N)  nth="$OPTARG";;
    ?)  usage;
    esac
done

if [ -z $staLIST ]
then
    echo "option -L must be specified!"
    usage
fi

# set parameters for RFmoho
CPTHOME=/usr/local/seis/data/cpt 

zmin=$(echo $zarg |cut -d"/" -f1)
zmax=$(echo $zarg |cut -d"/" -f2)
zinc=$(echo $zarg |cut -d"/" -f3)

rmin=$(echo $rarg |cut -d"/" -f1)
rmax=$(echo $rarg |cut -d"/" -f2)
rinc=$(echo $rarg |cut -d"/" -f3)

gridOUT=RFmoho.out

# set parameters for GMT
PS=RFmoho.ps
grdFile=RFmoho.grd
cptFile=$CPTHOME/RFmoho_lin.cpt
R=-R$zmin/$zmax/$rmin/$rmax
I=-I$zinc/$rinc
J=-JX5/5
B=-Ba10f5/a0.05f0.01WSne

# RFmoho
echo "RFmohoM -P $vp -Z $zmin/$zmax/$zinc -R$rmin/$rmax/$rinc -N$nth -L $staLIST -D $RFdir -w $w -t0.1 -n1 -s11 > $gridOUT"

#RFmohoM -P $vp -Z $zmin/$zmax/$zinc -R$rmin/$rmax/$rinc -N $nth -L $staLIST -D $RFdir -w $w -t0.1 -n1 -s11 -c > $gridOUT
RFmoho -P $vp -Z $zmin/$zmax/$zinc -R$rmin/$rmax/$rinc -N $nth -L $staLIST -D $RFdir -w $w -t0.1 -n1 -s11 > $gridOUT

# GMT
grep GRD $gridOUT | awk '{print $2, $3, $4*100.}' | xyz2grd -G$grdFile $I $R 
rrmax=$(grep HED $gridOUT |awk '{print $13}')
zzmax=$(grep HED $gridOUT |awk '{print $12}')

grdimage $grdFile $J -P -K  -V -C$cptFile $B -X1.5 -Y3.25  > $PS

# plot cross line
#psxy $R  $J -K -O -W5/255/255/255 -m <<END>>$PS
psxy $R  $J -K -O -W5/0/0/0 -m <<END>>$PS
<
$zzmax $rmin
$zzmax $rmax
<
$zmin $rrmax
$zmax $rrmax
END

psscale -C$cptFile -D2.5/-0.8/3.6/0.4h -O  -L  -B:." Beam Power ":  >> $PS
