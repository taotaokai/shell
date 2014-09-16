#!/bin/bash

function usage(){
    printf "Usage: %s: [-L SAC_List] [-D SAC_DIR] \n" $0
    echo "Default: $0 [-L-] -D."
	echo "evdp(km) and gcarc(degree) should be correctly set in SAC header!"
	echo "user0: rayp(s/km); user1: rayp(s/deg)"
    exit 1
}

# Default parameters

saclst=-
sacdir=.

# parse options
while getopts L:D:h name
do
    case $name in
    L)  saclst="$OPTARG";;
    D)  sacdir="$OPTARG";;
    [h,?])  usage; exit -1
    esac
done

#if [ -z "$saclst" ]
#then
	# read from stdin
#	saclst=-
    #echo "option -L must be specified!"
    #usage
#fi


#while read sacf
for sacf in $(cat $saclst)
do
	stninf=$(dumpSHD $sacdir/$sacf gcarc evdp)
	gcarc=$(echo $stninf | awk '{print $3}')
	evdp=$(echo $stninf | awk '{print $5}')

	angrayp=$(taup_time -h $evdp -deg $gcarc -ph P -rayp | sed -n "1p")
#	echo taup_time -h $evdp -deg $gcarc -ph P -rayp
#	echo "$angrayp * 180 / 6371 / 4 / a(1)" 
	horirayp=$(echo "scale=12;$angrayp * 180 / 6371 / 4 / a(1)" | bc -l)

	echo "$sacf gcarc $gcarc evdp $evdp angrayp $angrayp horirayp $horirayp"

sac<<EOF > /dev/null
r $sacdir/$sacf
*lh user0 user1
ch user0 $horirayp user1 $angrayp
lh user0 user1
wh
q
EOF

done
