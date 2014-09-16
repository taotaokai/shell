#!/bin/bash

#Convert FDSN catalog format into hypo format
#FDSN: 4723204|2014-05-27T07:52:23|44.6633|124.1042|14.4|US|NEIC PDE|NEIC COMCAT|1403854981000|MB|4.3|US|NORTHEASTERN CHINA
#hypo:         2014 06 17 00 33 39.0  44.6100  124.2600  15.0 3.7 ML|CENC 20140617_003339 eq 吉林前郭

CAT=$1

if [ -z "$CAT" ]; then
	echo "usage: $0 FDSN.cat"
	exit -1
fi

while IFS= read -r -u3 line
do

	char1=$(echo $line | sed "s/^[ \t]*//" | cut -c1)	
	if [ "$char1" == "#" ]; then
		echo $line
		continue
	fi

	other=$(echo $line | awk -F"|" '{print $13}')
	echo $line | awk -F"|" '{gsub(/[-:T]/," ",$2); print $2,$3,$4,$5,$11,$10}' | awk -v a="$other" '{printf "%04d %02d %02d %02d %02d %04.1f %8.4f %9.4f %5.1f %3.1f %3s|NEIC %04d%02d%02d_%02d%02d%02d %s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1,$2,$3,$4,$5,$6,a}'

done 3<$CAT

# END
