#!/bin/bash

#Synopsis
# RF_chkEVT.sh [dirnm]
#Description
#	check RF data under a given directory {dirnm}:
#	the paths of radial RFs are assumed to be {dirnm}/c_rfa/*.BHR
#	the direct P should be shifted to zero time
#	markers: bad(t5<0); needHPfilter(t6<0); reversedP(t7<0)
#History
#	[20140721]tao: created

if [ $# -ne 1 ]; then
	echo "USAGE: $0 [dirnm]"
	exit -1
fi

umask 002

wkdir=$(pwd)

dirnm=$1

echo "...processing $dirnm..."

RFdir=$wkdir/$dirnm/c_rfa/
if [ ! -d $RFdir ]; then 
	echo "$dirnm/c_rfa not exist"
	exit -1
fi	

cd $RFdir
nRF=$(ls *.BHR | wc -l | awk '{print $1}')
if [ $nRF -lt 1 ]; then
	echo "no BHR"
	exit -1
fi

sac<<END
r *.BHR
fileid l ur t l KSTNM KEVNM GCARC
qdp off
window x 0 1. y 0.2 1.
xlim -100 100
ppk p 10 bell off
wh
q
END

#	record selection
ls *.BHR | dumpSHD t5 | sed 's/\.BHR//' | awk '$3!=-12345 && $3<0 {print "'$evnm'",$1}' > bad.lst
ls *.BHR | dumpSHD t6 | sed 's/\.BHR//' | awk '$3!=-12345 && $3<0 {print "'$evnm'",$1}' > needHPfilter.lst
ls *.BHR | dumpSHD t7 | sed 's/\.BHR//' | awk '$3!=-12345 && $3<0 {print "'$evnm'",$1}' > reversedP.lst

#END
