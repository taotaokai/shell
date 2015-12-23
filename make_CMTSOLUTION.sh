#!/bin/bash

# create CMTSOLUTION file from GCMT ndk file for the given event gcmt-ID

#====== command line args
evid=${1:?evid not set}

outfile=${2:-${evid}.CMTSOLUTION}

ndk_file=${3:-GCMT_1976-NOW.ndk}

#====== 
tmp_ndk=$(mktemp)
grep -A 3 -B 1 "^$evid" $ndk_file > $tmp_ndk

sed '1!d;s/[\/:]/ /g;q' $tmp_ndk  > $outfile
sed '2!d;q' $tmp_ndk | awk '{printf "event name:     %s\n",$1}' >> $outfile
sed '3!d;q' $tmp_ndk | awk '{printf "time shift:     %8.4f\n",$2}' >> $outfile
sed '2!d;q' $tmp_ndk | awk '{printf "half duration:  %8.4f\n",$NF}' >> $outfile
sed '3!d;q' $tmp_ndk | awk '{printf "latitude:       %8.4f\n",$4}' >> $outfile
sed '3!d;q' $tmp_ndk | awk '{printf "longitude:      %8.4f\n",$6}' >> $outfile
sed '3!d;q' $tmp_ndk | awk '{printf "depth:          %8.4f\n",$8}' >> $outfile

sed '4!d;q' $tmp_ndk | awk '{printf "Mrr:            %8.4fe+%02d\n",$2,$1}' >> $outfile
sed '4!d;q' $tmp_ndk | awk '{printf "Mtt:            %8.4fe+%02d\n",$4,$1}' >> $outfile
sed '4!d;q' $tmp_ndk | awk '{printf "Mpp:            %8.4fe+%02d\n",$6,$1}' >> $outfile
sed '4!d;q' $tmp_ndk | awk '{printf "Mrt:            %8.4fe+%02d\n",$8,$1}' >> $outfile
sed '4!d;q' $tmp_ndk | awk '{printf "Mrp:            %8.4fe+%02d\n",$10,$1}' >> $outfile
sed '4!d;q' $tmp_ndk | awk '{printf "Mtp:            %8.4fe+%02d\n",$12,$1}' >> $outfile

rm $tmp_ndk
