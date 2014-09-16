#!/bin/bash

function usage(){
    printf "Usage: %s: [-C CAT_List] [-S STN_List] [-F from] [-T to] [-D DATA_DIR] [-O Out_DIR]\n" $0
    echo "Default: $0 -C CAT.lst -S- -B-50 -E500 -D. -O."
	echo "CAT_List: WEED file format (http://www.iris.edu/SeismiQuery/sq-events.htm)"
	echo "STN_List: netwk stnm stla stlo stel"
    echo "SEED data is stored in DATA_DIR/stnm/BH[E,N,Z]/(epochtime[ns]).seed"
	echo "DATA_DIR/stnm/SEED.lst should already exist"
    echo "SEED.lst: STN/BH[E,N,Z]/(epochtime[ns]).seed start-time stop-time"
    exit 1
}

# Default parameters

CATLST=CAT.lst
STNLST=-
FROM=-50
TO=500
OUTDIR=.
DATADIR=.

# parse options
while getopts C:S:F:T:D:O:h name
do
    case $name in
    C)  CATLST="$OPTARG";;
    S)  STNLST="$OPTARG";;
    F)  FROM="$OPTARG";;
    T)  TO="$OPTARG";;
    D)  DATADIR="$OPTARG";;
    O)  OUTDIR="$OPTARG";;
    [h,?])  usage; exit -1
    esac
done

if [ ! -f "$CATLST" ]
then
   echo "Error: cannot read CAT_LST!"
   usage
fi

# loop for each station

tmp_seed=$(mktemp -p.)
tmp_sac=$(mktemp -p.)
tmp_mac=$(mktemp -p.)

cat $STNLST | while read netwk stnm stla stlo stel tmp
do

    if [ ! -d $OUTDIR/$stnm ]
    then
        mkdir -p $OUTDIR/$stnm
    fi

	sed "s/,//g" $CATLST | while read tmp ymd mhs evla evlo evdp tmp tmp mag
	do
		# origin time (epoch time in Nano seconds)
		otime=$(date -u -d "$ymd $mhs" +%s%N)
		#ttime=$(aktimes $evla $evlo $evdp $stla $stlo | grep "[\t ]P[\t ]*$" | awk 'NR==1{print $1}')

		# travel-time&phase of the first arrival
		ttinfo=$(aktimes $evla $evlo $evdp $stla $stlo | head -n1)
		ttime=$(echo $ttinfo | awk '{print $1}')
		phase=$(echo $ttinfo | awk '{print $NF}')

		# first arrival time, begin&end cut time
		atime=$(echo "$otime+$ttime*10^9"|bc)
		btime=$(echo "$atime+$FROM*10^9"|bc)
		etime=$(echo "$atime+$TO*10^9"|bc)

		echo ===================
		echo Date: $(date)
		echo EVT: $ymd $mhs $mag
		echo STN: $netwk $stnm $stla $stlo $stel
        echo Phase: $phase $ttime $otime $atime $btime $etime

		# find seed files overlap time range [btime,etime]
		grep "$stnm" $DATADIR/$stnm/SEED.lst | awk '$3>=b && $2<=e {print dir"/"$0}' dir=$DATADIR e=$etime b=$btime > $tmp_seed
		if [ ! -s $tmp_seed ]
		then
			echo Warning: No file found!
			continue
		fi
		echo Files:
		cat $tmp_seed

		# convert seed to sac
		mseed2sac $(cat $tmp_seed | sort -k1 | awk '{print $1}') 2>$tmp_sac

		# time cut and write event info
		awk '{printf "r more %s\n",$NF}' $tmp_sac > $tmp_mac
		tmp=$(echo "$atime/10^9" | bc -l)
		gmt_atime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S.%N)
		tmp=$(echo "$otime/10^9" | bc -l)
		gmt_otime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S.%N)
		echo CUT: $gmt_otime, $gmt_atime
		echo "ch o GMT $gmt_otime" >> $tmp_mac
		echo "ch a GMT $gmt_atime" >> $tmp_mac
		mag1=$(echo $mag | grep MW | sed "s/.*MW[\t ]*\([^\t ]*\).*/\1/") # use MW magnitude
		if [ -z "$mag1" ] # if no MW found, use the last one
		then
			mag1=$(echo $mag | awk '{print $NF}')
		fi
		echo "ch knetwk $netwk kstnm $stnm stla $stla stlo $stlo stel $stel" >> $tmp_mac
		echo "ch evla $evla evlo $evlo evdp $evdp mag $mag1" >> $tmp_mac
		echo "wh" >> $tmp_mac
		echo "cut a $FROM $TO" >> $tmp_mac
		awk 'NR==1{printf "r %s\n",$NF} NR>=2{printf "r more %s\n",$NF}' $tmp_sac >> $tmp_mac
		echo "w over" >> $tmp_mac
sac <<EOF
m $tmp_mac
q
EOF
		# re-name sac and move into OUTDIR/stnm
		gmt_otime=$(date -u -d "$ymd $mhs" +%Y_%j_%H_%M_%S)
		awk '{printf "%s\n",$NF}' $tmp_sac |\
			awk -F"." '{printf "mv %s %s/%s/%s.%s.%s\n",$0,dir,stnm,stnm,o,$4}' o=$gmt_otime dir=$OUTDIR stnm=$stnm | bash
	done
done

# delete tmp files
rm $tmp_seed $tmp_sac $tmp_mac

# end of script
