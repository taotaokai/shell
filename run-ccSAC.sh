#!/bin/bash

#[20140802]tao: created

usage(){
cat <<EOF
Description:
	a utility script for ccSAC.

Usage: $0 -M<evtlst> -D<datlst> -O<outdir> -W<hd0>/<hd1>/<hd2>

	-M Specify the list of template wavform files (master file). 
		The list contains <path>/<evnm>.<netwk>.<stnm>.<comp>, e.g. <path>/20110510_023345.LX.SC12.BHZ
	-D Specify the list of data files.
		The list contains <path>/<evnm>.<netwk>.<stnm>.<comp>, e.g. <path>/20110510_023345.LX.SC12.BHZ
	-O Specify the directory for output, <outdir>/M<evnm>/D<datnm>/<netwk>.<stnm>.<comp>
	-W Specify the correlation time window on the template waveform.
		<hd0>: reference header name, e.g. o
		<hd1>: header name for the begin of cut window, e.g. t1
		<hd2>: header name for the end of cut window, e.g. t2 

Default: $0 -M evt.lst -D dat.lst -O.
EOF
	exit -1
}

# Default parameters
evtlst='evt.lst'
datlst='dat.lst'
outdir='.'
twin='o/t1/t2'

# parse options
while getopts M:D:O:W:h optname 
do
    case $optname in
    M)  evtlst="$OPTARG";;
    D)  datlst="$OPTARG";;
    O)  outdir="$OPTARG";;
    W)  twin="$OPTARG";;
    [h,?])  usage;;
    esac
done

# validate options
if [ ! -f "$evtlst" ] || [ ! -f "$datlst" ] || [ ! -d "$outdir" ]
then
	echo "Invalid option argument!"
	usage
fi

echo "#This run: $0 -M $evtlst -D $datlst -O $outdir -W $twin"

#
#tmplst=$(mktemp)

# loop each channel 
for evtfn in $(cat $evtlst)
do
	# get event name and channel name, e.g. 20110510_023345, LX.SC12.BHZ
	evnm=$(echo $evtfn | awk -F "/" '{print $NF}' | awk -F"." '{print $1}')
	chnm=$(echo $evtfn | awk -F "/" '{print $NF}' | awk -F"." '{printf "%s.%s.%s",$2,$3,$4}')

	# loop each data with the same channel name
	#grep "$chnm" $datlst > $tmplst
	for datfn in $(grep -F "$chnm" $datlst)
	do
		# get data name
		datnm=$(echo $datfn | awk -F "/" '{print $NF}' | awk -F"." '{print $1}')
		# output file name
		outfn_dir=$outdir/M$evnm/D$datnm
		outfn=$outfn_dir/$chnm

		# generate shell command
		echo
		echo "#M $evtfn"
		echo "#D $datfn"
		echo "#O $outfn"
		#echo cp -f $evtfn temp.sac
		#echo cp -f $datfn data.sac
		#echo chmod 755 temp.sac data.sac
		#printf "sac<<EOF\nr temp.sac data.sac;rmean;rtrend;taper;bp co 2 10 n 4 p 2; w over;q\nEOF\n"
		printf "if [ ! -d $outfn_dir ]; then\n\tmkdir -p $outfn_dir\nfi\n"
		#echo ccSAC -m temp.sac -d data.sac -o $outfn -w$twin
		echo ccSAC -m $evtfn -d $datfn -o $outfn -w$twin

	done
done

#rm $EVTLST $DATLST

# END
