#!/bin/bash

while read stnm
do

    mkdir $stnm/BHE $stnm/BHN $stnm/BHZ

    # make script for re-naming files
    echo "echo processing $stnm/BHE"
    grep "3.sac" $stnm.lst |\
        awk '{printf "mv %s %s/BHE/%s.sac\n",$1,stnm,$2}' stnm=$stnm
    echo "echo processing $stnm/BHN"
    grep "2.sac" $stnm.lst |\
        awk '{printf "mv %s %s/BHN/%s.sac\n",$1,stnm,$2}' stnm=$stnm
    echo "echo processing $stnm/BHZ"
    grep "1.sac" $stnm.lst |\
        awk '{printf "mv %s %s/BHZ/%s.sac\n",$1,stnm,$2}' stnm=$stnm

    # make data list of the re-named files
    grep 3.sac $stnm.lst |\
        awk '{printf "%s/BHE/%s.sac %s %s\n",stnm,$2,$2,$3}' stnm=$stnm > $stnm/SAC.lst
    grep 2.sac $stnm.lst |\
        awk '{printf "%s/BHN/%s.sac %s %s\n",stnm,$2,$2,$3}' stnm=$stnm >> $stnm/SAC.lst
    grep 1.sac $stnm.lst |\
        awk '{printf "%s/BHZ/%s.sac %s %s\n",stnm,$2,$2,$3}' stnm=$stnm >> $stnm/SAC.lst

done
