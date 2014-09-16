#!/bin/bash

tmarker=$1

if [ $# -eq 0 ];then
	tmarker='a'
fi

echo shift time to align "$tmarker" to zero

while read sacf
do

sac<<EOF
r $sacf
EVALUATE TO tmp &1,$tmarker * (-1)
ch allt %tmp
wh
q
EOF

done
