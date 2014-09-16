#!/bin/bash

# usage
function usage(){
	echo ==========================
	echo "Usage: $0 [-D DATADIR] [-I SACID] [-O OUTFILE]"
	echo --------------------------
	echo "DATADIR: where SAC files are stored"
	echo "SACID: find -name"SACID" will be used to locate all the sac files"
	echo "OUTFILE: output index file "
	echo ==========================
	echo "Default: $0 -D. -I\"*.sac\" -O SAC.index"
	echo ==========================
}

wkdir=$(pwd)

# Default parameters
DATDIR=.
SACID="*.sac"
OUTFN="SAC.index"

# parse options
while getopts D:I:O:h name
do
    case $name in
    D)  DATDIR="$OPTARG";;
    I)  SACID="$OPTARG";;
    O)  OUTFN="$OPTARG";;
    [h,?])  usage; exit 1;
    esac
done

# read sac file and output the start/end time in epoch nano-seconds
cat /dev/null > $wkdir/$OUTFN

cd $DATDIR
find -name "$SACID" | sed "s/\.\///" | sort |\
while read sacfile 
do
	gmt_ref=$(dumpSHD $sacfile nzyear nzjday nzhour nzmin nzsec nzmsec b e| awk '{print $3,$5,$7,$9,$11,$13,$15,$17}')
	str_ref=$(echo $gmt_ref | awk '{printf "%s/1/1 %s:%s:%s %d days",$1,$3,$4,$5+$6/1000,$2-1}')
	str_b=$(echo $gmt_ref | awk '{printf "%f",$7}')
	str_e=$(echo $gmt_ref | awk '{printf "%f",$8}')

	eptim_ref=$(date -u -d "$str_ref" +%s%N)
	eptim_b=$(echo "$eptim_ref+$str_b*10^9" | bc )
	eptim_e=$(echo "$eptim_ref+$str_e*10^9" | bc )
	echo $sacfile ${eptim_b%.*} ${eptim_e%.*} >> $wkdir/$OUTFN
done

# END
