#!/bin/bash

# plot and check SAC files 

if [ $# -ne 2 ] 
then
	echo "Two Arguments expected. Usage: chkSAC datlist tracenumber" 
	exit 1
fi

datlist=$1
goodlist=$1.good
badlist=$1.bad

list=$(cat $datlist | tr '\n' ' ')

echo $list
sac <<END > /dev/null
r $list
ch t5 1
wh
sort gcarc
fileid l ur t l KSTNM KEVNM GCARC
qdp off
ppk p $2 bell off
wh
q
END

if [ -f $goodlist ]
then
	cat /dev/null > $goodlist
fi

if [ -f $badlist ]
then
	cat /dev/null > $badlist
fi

for sac in $list
do
	dumpSHD $sac t5 | grep -v 12345 | awk '{if($3 > 0.) print $1}' >> $goodlist
done

cat $datlist $goodlist | awk '{a[$0]++} END{for ( x in a ) if ( a[x] == 1 ) print x}' >> $badlist

exit 1
