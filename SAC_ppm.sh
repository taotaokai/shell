#!/bin/bash

twin=$1

for sacf in $(ls *BHE | sed 's/\.BHE//')
do

sac<<EOF
cut $(echo $twin | sed "s/\// /")
r $sacf.BHN $sacf.BHE
bg sgf
ppm
sgftops f001.sgf $sacf.ps
q
EOF

done
