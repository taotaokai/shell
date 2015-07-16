#!/bin/bash

function usage(){
	echo ======================================
	printf "Usage: %s: [-C CUT_List] [-S STN_List] [-K CMP] [-D DATA_DIR] [-O Out_DIR]\n" $0
	echo --------------------------------------
	echo "CUT_List: evtnm start-time end-time"
	echo "For CUT_LIST, start-time and end-time should be in the format 1970-01-01T00:00:00 and timezone UTC+00:00 is assumed"
	echo "STN_List: netwk stnm stla stlo stel"
	echo "CMP: a string identifies components separated by slashes, such as .Z/.N/.E, each substring will be used with grep -F (treat as a fixed string!) to find out the corresponding files"
	echo "SAC data is stored in DATA_DIR/{stnm}/*.sac"
	echo "DATA_DIR/{stnm}/{stnm}.index should already exist"
	echo "The index file contains 3 columns: {stnm}/{sacfile} start-time end-time"
	echo "For the index file, start-time and end-time should be the epoch time (nanoseconds) equivalent to the output of \$(date -d "date string" +%s%N)"
	echo ======================================
	echo "Default: $0 -C CUT.lst -S- -K.Z/.N/.E -D. -O."
	echo ======================================
	exit 1
}

# Default parameters

CUTLST=CUT.lst
STNLST=-
CMPSTR=".E/.N/.Z"
OUTDIR=.
DATDIR=.

# parse options
while getopts C:S:K:D:O:h name
do
    case $name in
    C)  CUTLST="$OPTARG";;
    S)  STNLST="$OPTARG";;
    K)  CMPSTR="$OPTARG";;
    D)  DATDIR="$OPTARG";;
    O)  OUTDIR="$OPTARG";;
    [h,?])  usage; exit -1
    esac
done

if [ ! -f "$CUTLST" ]
then
   #echo "Error: option -L and -C must be specified!"
   echo "Error: $CUTLST does not exist!"
   usage
fi


# loop for each station
tmp_sac=$(mktemp -p.)
tmp_mac=$(mktemp -p.)
#
cat $STNLST | while read netwk stnm stla stlo stel tmp
do

    if [ ! -d $OUTDIR/$stnm ]
    then
        mkdir -p $OUTDIR/$stnm
    fi

	cat $CUTLST | while read starttime endtime
	do
			
		# begin/end cut time (nanoseconds)
		starttime=$(echo $starttime | sed "s/T/ /")
		endtime=$(echo $endtime | sed "s/T/ /")
		btime=$(date -u -d "$starttime" +%s%N )
		etime=$(date -u -d "$endtime" +%s%N )

		echo ===================
		echo Date: $(date)
		echo CUT: $starttime $endtime
		echo STN: $netwk $stnm $stla $stlo $stel

		# find sac files overlap time range [btime,etime]
		grep "$stnm" $DATDIR/$stnm/$stnm.index | awk '$3>=b && $2<=e {print dir"/"$0}' dir=$DATDIR e=$etime b=$btime | sort > $tmp_sac
		if [ ! -s $tmp_sac ]
		then
			echo Warning: No file found!
			continue
		fi
		echo Files:
		cat $tmp_sac

		# merge sac
		sacfn=$(mktemp -p. -u)
		cat /dev/null > $tmp_mac
		for cmpnm in $(echo $CMPSTR | sed "s/\// /g")
		do
			cat $tmp_sac | grep -F "$cmpnm" |\
       			     awk 'NR==1{printf "r %s\n",$1}NR>=2{printf "merge %s\n",$1}END{printf "w %s%s\n",sacfn,cmpnm}' sacfn="$sacfn" cmpnm="$cmpnm" >> $tmp_mac
		done

		# time cut and write event info
		tmp=$(echo "$btime/10^9" | bc -l)
		gmt_btime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S.%N)
		tmp=$(echo "$etime/10^9" | bc -l)
		gmt_etime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S.%N)
		echo "r $sacfn*" >> $tmp_mac
		echo "ch t1 gmt $gmt_btime" >> $tmp_mac
		echo "ch t2 gmt $gmt_etime" >> $tmp_mac
		echo "wh" >> $tmp_mac
		echo SAC: cut gmt $gmt_btime gmt $gmt_etime
		echo "cut t1 0 t2 0" >> $tmp_mac
		echo "r $sacfn*" >> $tmp_mac
		echo "ch knetwk $netwk kstnm $stnm stla $stla stlo $stlo stel $stel" >> $tmp_mac
		#echo "ch evla $evla evlo $evlo evdp $evdp mag $mag1" >> $tmp_mac
		echo "w over" >> $tmp_mac
sac <<EOF
m $tmp_mac
q
EOF

		# re-name sac and move to the output directory
		for cmpnm in $(echo $CMPSTR | sed "s/\// /g")
		do
			tmp=$(echo "$btime/10^9" | bc -l)
			gmt_btime=$(date -u -d @$tmp +%FT%T)
			echo mv ${sacfn}${cmpnm} $OUTDIR/$stnm/$stnm.${gmt_btime}${cmpnm}
			mv ${sacfn}${cmpnm} $OUTDIR/$stnm/$stnm.${gmt_btime}${cmpnm}
		done

	done
done

# delete tmp files
rm  $tmp_sac $tmp_mac

# end of script
