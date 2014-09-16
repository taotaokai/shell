#!/bin/bash

function usage(){
	echo ======================================
	echo "Description:"
	echo "	a wrap-up script of cmpAZ, which measures sensor orientation using P-wave particle motion polarity"	
	echo "	the paths of SAC data are assumed to be {DATDIR}/{STNM}/{SUBDIR}/*.{CMPNM}"
	echo "	the output cmpAZ.out is placed under {DATDIR}/{STNM}/"
	echo "	the direct P arrival should be shifted to zero time"
	echo ======================================
	echo "Usage: $0 [-D DATDIR] [-S STNLST] [-B SUBDIR] [-C CMPNM] [-W WINDOW]"
	echo 
	echo "-D DATDIR"
	echo "	the index file must be supplied, and the first column is the directory name {DIRNM}"
	echo "	this program will append ( | checked/skipped {date} ) to the end of each line after the corresponding directory is checked"
	echo "-S STNLST"
	echo "	a text list of the station directory name {STNM}; if not specified, read from stdin"
	echo "-B SUBDIR"
	echo "	sub-directory where the target sac files locate"
	echo "-C CMPNM"
	echo "	a string specifying 3-compoent names, e.g. BHN/BHE/BHZ (use slash to separate different components)"
	echo "-W WINDOW"
	echo "	specify the time window used for analyzing the polarity, for example -2/10 means 2 seconds before and 10 seconds after the zero time"
	echo ======================================
	echo "Default: $0 -D. -S- -B \"\" -C BHN/BHE/BHZ -W -2/10" 
	echo ======================================
	exit -1
}

# Default parameters
DATDIR='.'
STNLST='-'
SUBDIR=''
CMPNM='BHN/BHE/BHZ'
WINDOW='-2/10'

# parse options
while getopts D:S:B:C:W:h name
do
    case $name in
    D)  DATDIR="$OPTARG";;
    S)  STNLST="$OPTARG";;
    B)  SUBDIR="$OPTARG";;
    C)  CMPNM="$OPTARG";;
    W)  WINDOW="$OPTARG";;
    [h,?])  usage;;
    esac
done

echo "This run: $0 -D $DATDIR -S $STNLST -B \"$SUBDIR\" -C $CMPNM -W $WINDOW"

wkdir=$(pwd)

DATDIR=$(readlink -f $DATDIR)

#to find the sac files with the specified component names, use ls | grep 'BHN\|BHE\|BHZ'
CMPNM_grep=$(echo $CMPNM | awk -F"/" '{printf "%s\\|%s\\|%s",$1,$2,$3}')
#to delete the cmpnm, use sed 's/BHN//;s/BHE//;s/BHZ//'; this is to check if all 3 components exist
CMPNM_sed=$(echo $CMPNM | awk -F"/" '{printf "s/\\.%s//;s/\\.%s//;s/\\.%s//",$1,$2,$3}')

# for each station directory
cat $STNLST |\
while read stnm tmp
do

	cd $DATDIR/$stnm/$SUBDIR
	ls | grep "$CMPNM_grep" | sed "$CMPNM_sed" | awk '++a[$0]==3' > LIST #only use data with all 3 components

	cmpAZ -D . -L LIST -W $WINDOW -I 1 -C $CMPNM -v > cmpAZ.out	

done

# END
