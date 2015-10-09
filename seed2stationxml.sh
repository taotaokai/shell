#!/bin/bash

# IRIS seiscode:
# https://seiscode.iris.washington.edu/projects/stationxml-converter

jarfile=~/usr/seiscode/stationxml-converter-1.0.9.jar

for seed_file in "$@"
#while read seed_file
do

  xml_file=${seed_file%.*}.xml

  echo $seed_file to $xml_file

  java -jar $jarfile -p -so CDSN -x $seed_file -o $xml_file

done