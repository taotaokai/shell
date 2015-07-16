#!/bin/bash

# read sac file and output the start/end time in epoch nano-seconds

while read sacfile 
do
	gmt_ref=$(dumpSHD $sacfile nzyear nzjday nzhour nzmin nzsec nzmsec b e| awk '{print $3,$5,$7,$9,$11,$13,$15,$17}')
	str_ref=$(echo $gmt_ref | awk '{printf "%s/1/1 %s:%s:%s %d days",$1,$3,$4,$5+$6/1000,$2-1}')
	str_b=$(echo $gmt_ref | awk '{printf "%f",$7}')
	str_e=$(echo $gmt_ref | awk '{printf "%f",$8}')

	eptim_ref=$(date -u -d "$str_ref" +%s%N)
	eptim_b=$(echo "$eptim_ref+$str_b*10^9" | bc )
	eptim_e=$(echo "$eptim_ref+$str_e*10^9" | bc )

	#echo $sacfile $gmt_ref $str_ref $str_b $str_e eptimref $eptim_ref eptimb $eptim_b eptime $eptim_e
	echo $sacfile ${eptim_b%.*} ${eptim_e%.*}

done
