#!/bin/bash

#Convert epoch time(seconds since 19700101_000000) to 
#UTC time in yyyymmdd_HHMMSS[.NNN] format

while read epoch
do
	date -u -d @$epoch +%Y%m%d_%H%M%S.%N |\
    awk -F"." '{if($2!=0) printf "%s.%s\n",$1,substr($2,1,3); 
                else printf "%s\n",$1}'
done
