#!/bin/bash

#Format IRIS WEED file into hypo format defined by niu 
#USAGE: bash weed2hypo.sh < weedoriginmag-4.pl
#weed: ISC/ISC, 2010/05/31 19:51:55.0000, 11.16, 93.7, 127.9, 46, 703, MW,  6.5 
#hypo: yyyy mm dd hh mm ss.sss evla evlo evdp mag type|origin evnm
#[2014-06-24]tao: created

while read line
do
	echo $line | sed "s/^ \t*//" |\
	awk -F"," '{a=substr($1,1,4);printf "%s %8.4f %9.4f %5.1f %3.1f %5s|%4s\n",$2,$3,$4,$5,$9,$8,a}' |\
	sed "s/\// /g;s/:/ /g" |\
	awk '{printf "%s  %s%s%s_%s%s%02.0f\n",$0,$1,$2,$3,$4,$5,$6}' |\
	sed "s/^ *//"
done
