#!/bin/bash

#Convert string of utc time (yyyymmdd_HHMMSS[.SSS]) to epoch time (second)

while read utc 
do
	ymd=${utc:0:4}-${utc:4:2}-${utc:6:2}
	hms=${utc:9:2}:${utc:11:2}:${utc:13}
	date -u -d "${ymd} ${hms}" +%s%N | awk '{print $1/10^9}'
done
