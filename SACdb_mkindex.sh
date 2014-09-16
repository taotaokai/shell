#!/bin/bash

# usage
function usage(){
	echo ==========================
	printf "Usage: %s: [-S STN_List] [-D DATA_DIR] [-I SACID] [-O Out_DIR]\n" $0
	echo --------------------------
	echo "STN_List: stnm"
	echo "DATA_DIR: where SAC files are stored"
	echo "SACID: find -name"SACID" will be used to locate all the sac files"
	echo "OUT_DIR: the index file stnm.index will be written into this dicretory"
	echo ==========================
	echo "Default: $0 -S- -I\"*.sac\" -D. -O."
	echo ==========================
}


wkdir=$(pwd)

# Default parameters
STNLST=-
DATDIR=.
SACID="*.sac"
OUTDIR=.

# parse options
while getopts S:D:I:O:h name
do
    case $name in
    S)  STNLST="$OPTARG";;
    D)  DATDIR="$OPTARG";;
    I)  SACID="$OPTARG";;
    O)  OUTDIR="$OPTARG";;
    [h,?])  usage; exit 1;
    esac
done

cat $STNLST |\
while read stnm
do
	echo =================
	echo Processing $stnm

	cd $DATDIR
	# read sac file and output the start/end time in epoch nano-seconds
	cat /dev/null > $wkdir/$OUTDIR/$stnm.index
  find $stnm -name "$SACID" | sort |\
	while read sacfile 
	do
		gmt_ref=$(dumpSHD $sacfile nzyear nzjday nzhour nzmin nzsec nzmsec b e| awk '{print $3,$5,$7,$9,$11,$13,$15,$17}')
		str_ref=$(echo $gmt_ref | awk '{printf "%s/1/1 %s:%s:%s %d days",$1,$3,$4,$5+$6/1000,$2-1}')
		str_b=$(echo $gmt_ref | awk '{printf "%f",$7}')
		str_e=$(echo $gmt_ref | awk '{printf "%f",$8}')

		eptim_ref=$(date -u -d "$str_ref" +%s%N)
		eptim_b=$(echo "$eptim_ref+$str_b*10^9" | bc )
		eptim_e=$(echo "$eptim_ref+$str_e*10^9" | bc )
		echo $sacfile ${eptim_b%.*} ${eptim_e%.*} >> $wkdir/$OUTDIR/$stnm.index
	done
	cd $wkdir
done
