#!/bin/bash

chnm=BH
if [ -n "$1" ]
then
  chnm=$1
fi

while read sacf
do

echo $sacf

sac<<EOF > /dev/null
r $sacf.${chnm}E
ch cmpaz 90 cmpinc 90
wh

r $sacf.${chnm}N
ch cmpaz 0 cmpinc 90
wh

r $sacf.${chnm}E $sacf.${chnm}N 
rotate 
w $sacf.${chnm}R $sacf.${chnm}T

r $sacf.${chnm}R; ch kcmpnm ${chnm}R; wh
r $sacf.${chnm}T; ch kcmpnm ${chnm}T; wh
q
EOF

done
