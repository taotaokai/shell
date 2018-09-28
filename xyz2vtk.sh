#!/bin/bash

input_xyz=${1:?[arg]need xyz list}

npoint=$(awk 'NF==3&&$1!~/#/' $input_xyz | wc -l)

cat<<EOF
# vtk DataFile Version 2.0
Source and Receiver VTK file
ASCII
DATASET POLYDATA
POINTS   $npoint   float
EOF

awk 'NF==3&&$1!~/#/{print $0}; END{printf"\n"}' $input_xyz
