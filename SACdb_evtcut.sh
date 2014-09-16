#!/bin/bash

#History:
#[20140624]tao: modified from SACdb_wincut.sh 
#[20140721]tao: change order "stla stlo" -> "stlo stla" 
#[20140721]tao: change default -K option to BHE/BHN/BHZ 
#[20140726]tao: used file descriptor to safely read lines from "while read -r"
#[20140726]tao: set origin time to zero
#[20140726]tao: remove the default behavior to read from stdin if -S not specified 

usage(){
	echo "Usage: $0 -E<hypo> -S<sta> -K<comp> -I<index> -W<from>/<to> -D<dataDIR> -O<outDIR>"
	echo
	echo "	-E Specify list of hypocenters. Each line containes: yyyy mm dd HH MM SS.S evla evlo evdp mag magType evnm"
	echo "	-S Specify list of stations. Each line contains: netwk stnm stlo stla stel"
	echo "	-K Specify component namess. Multiple components should be separated by slashes, such as BHE/BHN/BHZ, each substring will be used with grep -F (treat as a fixed string!) to find out the corresponding files of the same component for the merging process;"
	echo "	-I Specify index file. The index file contains (sacfn start-time end-time) on each line for every sac files in the data directory."
	echo "	   start-time and end-time are epoch time (nanoseconds) as the output of \$(date -d "date string" +%s%N)"
	echo "	-W Specify the cut time window. Time window (second) is relative to the first arrival."
	echo "	-D Specify the data directory."
	echo "	-O Specify the output directory."
	echo
	echo "Path of the output files: outDIR/stnm/evnm.netwk.stnm.comp"
	echo 
	echo "Default: $0 -E hypo -S sta -K BHZ/BHN/BHE -I SAC.index -W -100/500 -D. -O."
	exit -1
}

# Default parameters

EVTLST=hypo
STNLST=sta
CMPSTR=BHE/BHN/BHZ
INDEX=SAC.index
CUTWIN=-100/500
DATDIR=.
OUTDIR=.

# parse options
while getopts E:S:K:I:W:D:O:h name
do
    case $name in
    E)  EVTLST="$OPTARG";;
    S)  STNLST="$OPTARG";;
    K)  CMPSTR="$OPTARG";;
    I)  INDEX="$OPTARG";;
    W)  CUTWIN="$OPTARG";;
    D)  DATDIR="$OPTARG";;
    O)  OUTDIR="$OPTARG";;
    [h,?])  usage; exit -1
    esac
done

if [ ! -f "$EVTLST" ]
then
   echo "Error: hypo file <$EVTLST> does not exist!"
   usage
fi

if [ ! -f "$INDEX" ]
then
	echo "Warning: index file <$INDEX> does not exist, use $DATDIR/SAC.index"
	INDEX=$DATDIR/SAC.index
fi

if [ ! -f "$INDEX" ]
then
	echo "Error: cannot find $DATDIR/SAC.index!"
	usage
fi

echo "This run: $0 -E $EVTLST -K $CMPSTR -I $INDEX -W $CUTWIN -D $DATDIR -O $OUTDIR"

# time window
FROM=$(echo $CUTWIN | awk -F"/" '{print $1}')
TO=$(echo $CUTWIN | awk -F"/" '{print $2}')

# loop for each station
#tmp_sac=$(mktemp -p.)
#tmp_sacm=$(mktemp -p.)
tmp_lst=$(mktemp -p. tmp.lst.XXX)
tmp_sacm=$(mktemp -p. tmp.sacm.XXX)
tmp_sta=$(mktemp -p. tmp.sta.XXX)
tmp_evt=$(mktemp -p. tmp.cut.XXX)

# sed: remove blank and commented lines
sed '/^\s*$/d;/^#/d' $EVTLST > $tmp_evt
sed '/^\s*$/d;/^#/d' $STNLST > $tmp_sta

while read -r -u3 netwk stnm stlo stla stel tmp
do

	if [ ! -d $OUTDIR/$stnm ]
	then
	    mkdir -p $OUTDIR/$stnm
	fi

	while read -r -u4 year month day hour min sec evla evlo evdp mag magtype evnm junk 
	do
			
		# origin time
		odate="$year-$month-$day $hour:$min:$sec"
		otime=$(date -u -d "$odate" +%s%N)

		# travel-time of the first arrival
		ttinfo=$(aktimes $evla $evlo $evdp $stla $stlo ak135 | head -n1)
		ttime=$(echo $ttinfo | awk '{print $1}')
		phase=$(echo $ttinfo | awk '{print $NF}')

		# first arrival time, begin&end cut time
		atime=$(echo "$otime+$ttime*10^9"|bc)
		btime=$(echo "$atime+$FROM*10^9"|bc)
		etime=$(echo "$atime+$TO*10^9"|bc)

		echo ===================
		echo Date: $(date)
		echo EVT: $odate $evla $evlo $evdp $mag $magtype $evnm
		echo STN: $stnm $phase $ttime 
		echo CUT: o $otime a $atime b $btime e $etime

		# find sac files overlap time range [btime,etime] in the INDEX file
		grep "$stnm" $INDEX | awk '$3>=b && $2<=e {print dir"/"$0}' dir=$DATDIR e=$etime b=$btime | sort > $tmp_lst
		if [ ! -s $tmp_lst ]
		then
			echo Warning: No file found!
			continue
		fi
		echo Files:
		cat $tmp_lst

		# generate sac macro
		cat /dev/null > $tmp_sacm # the content from the last loop must be cleaned

		# merge sac
		sacfn=$(mktemp -p. -u)
		for cmpnm in $(echo $CMPSTR | sed "s/\// /g")
		do
			cat $tmp_lst | grep -F "$cmpnm" |\
 			awk 'NR==1{printf "r %s\n",$1}NR>=2{printf "merge %s\n",$1}END{if(NR>=1) printf "ch kcmpnm %s;w %s.%s\n",cmpnm,sacfn,cmpnm}' sacfn="$sacfn" cmpnm="$cmpnm" >> $tmp_sacm
		done

		# time cut and write event info
		tmp=$(echo "$atime/10^9" | bc -l)
		gmt_atime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S" "%N | awk '{printf "%s %s %s %s %s %3d",$1,$2,$3,$4,$5,$6/10^6}')
		tmp=$(echo "$otime/10^9" | bc -l)
		gmt_otime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S" "%N | awk '{printf "%s %s %s %s %s %3d",$1,$2,$3,$4,$5,$6/10^6}')

		echo SAC: cut $gmt_atime From $FROM To $TO

		echo "r $sacfn*" >> $tmp_sacm
		echo "ch o gmt $gmt_otime" >> $tmp_sacm
		echo "ch a gmt $gmt_atime" >> $tmp_sacm
		echo "wh" >> $tmp_sacm

		echo "cut a $FROM a $TO" >> $tmp_sacm
		echo "r $sacfn*" >> $tmp_sacm
		echo "ch kevnm $evnm knetwk $netwk kstnm $stnm stla $stla stlo $stlo stel $stel" >> $tmp_sacm
		echo "ch evla $evla evlo $evlo evdp $evdp mag $mag" >> $tmp_sacm
		echo "w over" >> $tmp_sacm

		#set origin time to zero 
		echo "cut off" >> $tmp_sacm
		for cmpnm in $(echo $CMPSTR | sed "s/\// /g")
		do
			echo "r ${sacfn}.${cmpnm}" >> $tmp_sacm
			echo "eval to v -1 * &1,o" >> $tmp_sacm
			echo "ch allt %v" >> $tmp_sacm
			echo "wh" >> $tmp_sacm
		done

sac <<EOF
m $tmp_sacm
q
EOF

		# re-name sac and move to the output directory
		for cmpnm in $(echo $CMPSTR | sed "s/\// /g")
		do
			echo mv ${sacfn}.${cmpnm} $OUTDIR/$stnm/${evnm}.${netwk}.${stnm}.${cmpnm}
			mv ${sacfn}.${cmpnm} $OUTDIR/$stnm/${evnm}.${netwk}.${stnm}.${cmpnm}
		done

	done 4<$tmp_evt
done 3<$tmp_sta

# delete tmp files
#
rm  $tmp_lst $tmp_sacm $tmp_evt $tmp_sta

# END of script
