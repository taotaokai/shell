#!/bin/bash

wkdir=$(pwd)

while read datdir stnm
do
	if [ ! -d $stnm ]
	then
		mkdir -p $wkdir/$stnm
	fi

	cd $wkdir/$stnm
    pwd
    ckreftek $datdir >> ckreftk.log

done
