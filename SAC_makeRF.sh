#!/bin/bash

while read evnm
do

sac<<EOF
r $evnm.BHR $evnm.BHZ
rmean;rtrend;
w BHR BHZ
cut -50 100
r BHR BHZ
w over
q
EOF

FDdeconv -S BHZ -F BHR -f-50 -t100 -s50 -k 0.01 -G -a2.5 -o $evnm.RZ 
sac<<EOF
r $evnm.RZ
ch b -50 
wh
q
EOF

done
