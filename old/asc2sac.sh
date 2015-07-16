#!/bin/bash

#for file in "$@"
while read file
do
	echo $file
	npts=$(wc -l $file|awk '{print $1}')
	asc2sac $file $npts $file.sac
done
