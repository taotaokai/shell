#!/bin/bash

#Change time zone in CENC catalog (GMT+8) to UTC(GMT+0)

#2014-06-17  08:33:39.0  44.61  124.26  15  ML  3.7  eq  吉林前郭

CAT=$1

while IFS= read -r -u3 line 
do
	
	char1=$(echo $line | sed "s/^[ \t]*//" | cut -c1)	
	if [ "$char1" == "#" ]; then
		echo $line
		continue
	fi

	gmtstr=$(echo $line | awk '{print $1,$2}')
	other=$(echo $line  | awk '{printf "%s %s %2d %s %3.1f %s %s",$3,$4,$5,$6,$7,$8,$9}')

	gmtstr2=$(date -u -d "$gmtstr +8" +"%Y-%m-%d %H:%M %S.%N" | awk '{printf "%s %s:%04.1f",$1,$2,$3}')

	echo "$gmtstr2 $other"

done 3<$CAT
