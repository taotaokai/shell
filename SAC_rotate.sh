#!/bin/bash

while read sacf
do

sac<<EOF
r $sacf.BHE
ch cmpaz 90 cmpinc 90
wh

r $sacf.BHN
ch cmpaz 0 cmpinc 90
wh

r $sacf.BHE $sacf.BHN 
rmean;rtrend;taper
hp co 0.01 n 4 p 2
rotate 
w $sacf.BHR $sacf.BHT

r $sacf.BHR; ch kcmpnm BHR; wh
r $sacf.BHT; ch kcmpnm BHT; wh
q
EOF

done
