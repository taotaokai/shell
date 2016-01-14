#!/bin/bash

selfdoc(){
cat <<EOF
NAME
  checkpair

DESCRIPTION
  check the content from stdin to see if paired

SYNOPSIS
  checkpair < sac_files -n 2

COMMENTS
  matched pairs will be printed out

PARAMETERS
  -r  if issued, only the unmatch pair will be printed out together
      with the number of occurrence
  -n 2  specify the number in the pair

EOF
}

# check if stdin(0) is opened
if [ -t 0 ]
then
  selfdoc 
fi

# default
n=2
flag_reverse=n

# parse options 
while getopts n:rh name
do
  case $name in
  n) n="$OPTARG";;
  r) flag_reverse=y;;
  [h,?]) selfdoc; exit -1
  esac
done

# do the real job 
if [ $flag_reverse = "y" ]
then
  cat - |\
  awk '{++a[$0]}; END{for(x in a) if(a[x]!=n) print x, a[x]}' n=$n
else
  cat - |\
  awk '{++a[$0]}; END{for(x in a) if(a[x]==n) print x}' n=$n
fi

#END
