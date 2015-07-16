#!/bin/bash

while read stnm
do
    find $stnm -name "*.sac" | sort > $stnm.sac.lst

    cat $stnm.sac.lst | bash SAC_gettime.bash > $stnm.lst

done
