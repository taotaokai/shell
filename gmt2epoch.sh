#!/bin/bash

#Convert string of gmt time (yyyymmdd_HHMMSS[.SSS]) to epoch time (second)

while read gmt
do
	ymd=${gmt:0:4}-${gmt:4:2}-${gmt:6:2}
	hms=${gmt:9:2}:${gmt:11:2}:${gmt:13}
	date -u -d "${ymd} ${hms}" +%s%N | awk '{print $1/10^9}'
done
