#!/bin/bash

while read sacfn 
do

	gmt_ref=$(dumpSHD $sacfn nzyear nzjday nzhour nzmin nzsec nzmsec b e| awk '{print $3,$5,$7,$9,$11,$13,$15,$17}')
	str_ref=$(echo $gmt_ref | awk '{printf "%s/1/1 %s:%s:%s %d days",$1,$3,$4,$5+$6/1000,$2-1}')
	str_b=$(echo $gmt_ref | awk '{printf "%f",$7}')
	str_e=$(echo $gmt_ref | awk '{printf "%f",$8}')

	#eptim_ref=$(date -u -d "$str_ref" +%s%N) #nano second
	eptim_ref=$(date -u -d "$str_ref" +%s) #second
	eptim_b=$(echo "$eptim_ref+$str_b" | bc )
	eptim_e=$(echo "$eptim_ref+$str_e" | bc )
	echo $sacfn b ${eptim_b%.*} e ${eptim_e%.*}

done
