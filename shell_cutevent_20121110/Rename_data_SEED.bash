#!/bin/bash

while read stnm
do

    mkdir $stnm/BHE $stnm/BHN $stnm/BHZ

    # make script for re-naming files
    echo "echo processing $stnm/BHE"
    grep BHE $stnm.lst |\
        awk '{printf "mv %s %s/BHE/%s.seed\n",$1,stnm,$2}' stnm=$stnm
    echo "echo processing $stnm/BHN"
    grep BHN $stnm.lst |\
        awk '{printf "mv %s %s/BHN/%s.seed\n",$1,stnm,$2}' stnm=$stnm
    echo "echo processing $stnm/BHZ"
    grep BHZ $stnm.lst |\
        awk '{printf "mv %s %s/BHZ/%s.seed\n",$1,stnm,$2}' stnm=$stnm

    # make data list of the re-named files
    grep BHE $stnm.lst |\
        awk '{printf "%s/BHE/%s.seed %s %s\n",stnm,$2,$2,$3}' stnm=$stnm > $stnm.inf
    grep BHN $stnm.lst |\
        awk '{printf "%s/BHN/%s.seed %s %s\n",stnm,$2,$2,$3}' stnm=$stnm >> $stnm.inf
    grep BHZ $stnm.lst |\
        awk '{printf "%s/BHZ/%s.seed %s %s\n",stnm,$2,$2,$3}' stnm=$stnm >> $stnm.inf

done
