#!/bin/bash

#Convert CENC catalog format into hypo format
#CENC: 2014-06-17 00:33:39.0 44.61 124.26 15 ML 3.7 eq 吉林前郭
#hypo: 2014 06 17 00 33 39.0  44.6100  124.2600  15.0 3.7 ML|CENC 20140617_003339 eq 吉林前郭

CAT=$1

if [ -z "$CAT" ]; then
	echo "usage: $0 CENC.cat"
	exit -1
fi

while IFS= read -r -u3 line
do

	char1=$(echo $line | sed "s/^[ \t]*//" | cut -c1)	
	if [ "$char1" == "#" ]; then
		echo $line
		continue
	fi

	other=$(echo $line | awk '{print $8,$9}')
	echo $line | sed "s/-/ /g;s/:/ /g" | awk -v a="$other" '{printf "%04d %02d %02d %02d %02d %04.1f %8.4f %9.4f %5.1f %3.1f %3s|CENC %04d%02d%02d_%02d%02d%02d %s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$10,$1,$2,$3,$4,$5,$6,a}'

done 3<$CAT

# END
