#!/bin/bash

selfdoc(){
cat <<EOF
NAME
  SAC_shift2zero.sh

DESCRIPTION
  shift sac time markers to make one header value zero.

SYNOPSIS
  SAC_shift2zero.sh < sac-files <tmarker>

PARAMETERS
  <tmarker>  specify the sac header name

EOF
  exit -1
}

# Default parameters
HEAD=${1:?[arg] need header name (e.g. o, a, t1, etc)}

# parse options
#while getopts t:h name
#do
#  case $name in
#  t)  HEAD="$OPTARG";;
#  [h,?])  selfdoc;;
#  esac
#done

echo shift sac to align header \""$HEAD"\" to zero

# loop each sac file 
while read sacf
do

sac<<EOF > /dev/null
r $sacf
EVALUATE TO tmp &1,$HEAD * (-1)
ch allt %tmp
wh
q
EOF

done

#END
