#!/bin/bash

#USAGE: gmt_chTZ.sh "-8" < time.lst (convert Beijing time to UTC)

diffT=$1

while read gmtstr
do

epoch1=$(date -u -d "$gmtstr" +%s.%N)
epoch2=$(echo $epoch1 + ${diffT}*3600 | bc)

gmtstr2=$(date -u -d @$epoch2 +"%Y%m%d_%H%M%S.%N")

echo $gmtstr2

done
