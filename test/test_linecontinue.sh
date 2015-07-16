#!/bin/bash

# Test various line continuation

text="some text"
var=end

echo $text |\
awk 'BEGIN{print "begin"}\
     {print $0}\
     END{print x}' \
     x=$var

echo "this is\
 a test"

cat <<EOF
this is cat
another line
EOF
