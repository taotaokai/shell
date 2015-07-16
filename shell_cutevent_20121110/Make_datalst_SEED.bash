#!/bin/bash

while read stnm
do
    find $stnm -name "*.m" | sort > $stnm.mseed.lst

    cat $stnm.mseed.lst | bash Seed_gettime.bash > $stnm.lst

done
