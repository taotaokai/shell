#!/bin/bash

while read year month day
do
	# convert from y,m,d to y,jday
	comm_bc1="a=(14-$month)/12;y=$year+4800-a;m=$month+12*a-3;$day+(153*m+2)/5+365*y+y/4-y/100+y/400-32045"
	comm_bc0="a=(14-1)/12;y=$year+4800-a;m=1+12*a-3; (153*m+2)/5+365*y+y/4-y/100+y/400-32045"
	jday1=$(echo "$comm_bc1" | bc)
	jday0=$(echo "$comm_bc0" | bc)
	jday=$(echo "$jday1-$jday0" | bc)

	echo $year $jday

done
