#!/bin/bash

# read seed file and output the start/end time in epoch nano-seconds

while read seedfile
do

	gmt_be=$(mseedhdr $seedfile -s 2>/dev/null | tail -n1 |  awk '{print $3,$5}')

    # skip broken or empty mseed file
    if [ -z "$gmt_be" ]
    then
        continue
    fi
 
	str_b=$(echo $gmt_be | awk '{print $1}' | awk -F":" '{printf "%s/1/1 %s:%s:%s %d days",$1,$3,$4,$5,$2-1}')
	str_e=$(echo $gmt_be | awk '{print $2}' | awk -F":" '{printf "%s/1/1 %s:%s:%s %d days",$1,$3,$4,$5,$2-1}')

	eptim_b=$(date -u -d "$str_b" +%s%N)
	eptim_e=$(date -u -d "$str_e" +%s%N)

	#echo $seedfile $gmt_be $str_b $str_e $eptim_b $eptim_e
	echo $seedfile $eptim_b $eptim_e

done
